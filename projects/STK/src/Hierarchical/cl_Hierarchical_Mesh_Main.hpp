/*
 * cl_Hierarchical_Mesh_Main_Main.hpp
 *
 *  Created on: Dec 7, 2017
 *      Author: gleim
 */

#ifndef SRC_MESH_CL_Hierarchical_Mesh_Main_MAIN_HPP_
#define SRC_MESH_CL_Hierarchical_Mesh_Main_MAIN_HPP_
// MORIS library header files.
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
#include "fn_linsolve.hpp"
#include "fn_save_vector_to_binary_file.hpp" // LNA/src
#include "cl_Solver_Input.hpp" // DLA/src
// #include "cl_Bitset.hpp" // CON/src
#include "cl_Map.hpp" // CON/src
#include "chronos.hpp"
#include "cl_Communication_Tools.hpp" // COM/src
#include "cl_Bspline.hpp" // MOD/src
#include "cl_BoostBitset.hpp" // CON/src
#include "cl_Mesh_Enums.hpp" // MTK/src

#include "cl_Hierarchical_Mesh_Element.hpp"
#include "cl_Hierarchical_Mesh_Basis.hpp"
#include "cl_Hierarchical_Mesh_Edge.hpp"
#include "cl_Hierarchical_Mesh_Face.hpp"
#include "cl_Hierarchical_Mesh_TMatrix.hpp"
#include "cl_Hierarchical_Mesh_Refinement.hpp"
#include "cl_Hierarchical_Mesh_Sets.hpp"
#include "cl_Hierarchical_Mesh_MPI.hpp"
#include "cl_Hierarchical_Mesh_MTK.hpp"
#include "cl_Hierarchical_Mesh_Output.hpp"
#include "cl_Hierarchical_Mesh_Input.hpp"
#include "cl_Hierarchical_Mesh_Filter.hpp"

#include "cl_Lagrange_Basis.hpp"
#include "cl_Lagrange_Element.hpp"
#include "cl_Lagrange_Edge.hpp"
#include "cl_Lagrange_Face.hpp"
#include "cl_Lagrange_Filter.hpp"

#include "cl_Base_Mesh_Element.hpp"
#include "cl_Base_Mesh_Edge.hpp"

#include "cl_Solver_Factory.hpp" // DLA/src

namespace moris
{
    class Hierarchical_Mesh_Solver_Input;

// -----------------------------------------------------------------------------

    // User defined settings
    struct HMR_Settings
    {
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // General Settings
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        // flag telling if timing information should be printed
        bool TimingTest = false;

        // flag to check if connectivity data are to be interpreted and written
        // as binary
        bool InputOutputFilesForMorisAreAscii = false;

        // flag to check if connectivity data are to be interpreted and written
        // as binary
        bool OutputFilesForFemdocAreBinary = false;

        // probably deprecated function. This switch is always false.
        bool DesignVariablesOnInitialMesh = false;

        // this function is currently broke. FIXME: fix this!
        bool RefineBaseOnExoFile = false;

        // flag to tell if design variables are to be projected on output mesh
        bool ProjectionFlag = false;

        // coordinate of first node in user specified mesh (without invisible elements)
        Mat<real> PointOfOrigin = {{0},{0},{0}};

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // Refinement Settings
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        // Flag telling if B-Splines are to be truncated (-->partition of unity)
        bool TruncatedBsplines = false;
        // deprecated feature.
        bool PerformNormalization = false;

        // switch telling if current step is an initialization step (false)
        // or refinement step (true)
        bool Refinement = false;

        // all elements are to be refined to at least this level
        uint MinimumRefinement = 0;

        // minimum number of neighbor elements to be considered in buffer layer
        uint BufferElements = 0;

        // maximum refinement level for element field (non-SDF)
        uint MaxDesignLevelOfRefinement = 3;

        // maximum refinement level for surface field (non-SDF)
        uint MaxLevelSetLevelOfRefinement = 3;

        // lower bound for refining against density field
        real ElementLowerBound = 0.1;
        real ElementIncrementBound = 0.0; // deprecated parameter
        real ElementUpperBound = 1.0;

        real NodalLowerBound = 0.0;
        real NodalUpperBound = 0.0;

        bool AdaptiveBufferLayer = false;
        bool Staircasebuffer = true;

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // Settings for SDF
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        // flag if SDF is to be used
        bool SDFField = false;

        // calculate SDF for surfaces for the next n Neighbors
        uint SDFCandidateSearchDepth = 0;

        // cell containing paths to wavefront object files for SDF
        Cell< std::string > SDFObjectFiles;

        // maximum refinement level for SDF surface
        uint MaxSurfaceRefinement = 0;

        // maximum refinement level for SDF volume
        uint MaxVolumeRefinement = 0;

        // scaling parameter for level set creation out of SDF
        Mat< real > SDFAlpha = { {1}, {1}, {1}, {1}, {1}, {1} };

        // offset parameter for level set creation out of SDF
        Mat< real > SDFBeta  = { {0}, {0}, {0}, {0}, {0}, {0} };

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // Settings for Level Set Creation
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        // switch used in find_elements_for_refinement (FIXME: deprecated?)
        bool LvlSetMethod = false;

        // exponent for SIMP density field
        real SimpExp = 1.0;

        // epsilon environment near surface defining if a point is to be
        // associated with SDF surface or volume
        real ElementEdgeWidth        = 1;
        real InitialElementEdgeWidth = -1; // negative number indicates that ElementEdgeWidth should be used

        // relative scaling factor for level set
        // ( absolute factor is LSscale*ElementEdgeWidth )
        real LSscale  = 5.0;

        // level set threshold
        real LSthresh = 0.5;

        // special parameter for level set creation
        real DensityShift = 0.01;

        // special parameter for level set creation
        real Prjbeta = 0.001;

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // Settings for Feature size control
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        uint FeatureResolution = 3;
        real FeatureLowerBound = 0.5; // FIXME: is this parameter actually used?

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // Filter
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        real FilterRadius      = 0.0;
        bool FilterLevelWeight = false;

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // Symmetry
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        // flag telling if symmetry is to be enforced
        bool UseSymmetry = false;

        // first value: 0: x-axis, 1: y-axis, 2: z-axis
        uint  SymmetryPlane =  2 ;

        // radius for holes to be generated. Off if radius equal 0
        real HoleRadius = 0;

        real HoleExponent = 2;

    };

// -----------------------------------------------------------------------------

    struct HMR_Element_Data
    {
        Mat<uint> DeactivateElement;    // list of elements which are to be refined in the next step

        BoostBitset ElementActive;
        BoostBitset ElementRefined;

        BoostBitset ElementActiveDummy;    // temporary variable needed for MTK output
        BoostBitset ElementActiveDesign;
        BoostBitset ElementRefinedDesign;
        BoostBitset ElementActiveLastStep;
        BoostBitset ElementActiveDesignLastStep;

        Mat<uint> ElementListActiveDesign;

        Mat<uint> FeTopo;
        Mat<uint> FeTopoMTK;
        Mat<real> FeCoord;

        Mat< real > TMatrix;
        Mat< uint > IdField;
        Mat< real > TMatrixDesign;
        Mat< uint > IdFieldDesign;

        Cell<Mat<real>> TMatrixParentChildRelationFEM;
        Cell<Mat<real>> TMatrixParentChildRelationDesign;

        // this one can be calculated in create_mesh, but is never used
        //Cell<Mat<real>> TMatrixParentChildRelationDesignProjection;

        Mat<real> TProjectionMatrixBasis;
        Mat<real> TProjectionMatrixDesignBasis;

        // undocumented items below

        // FIXME: this one is used in create_mesh to make sure that the max number of possible IDs does not
        // exceed UINT_MAX. But the data is not consistent. This needs to be done in another way.
        Mat<uint> MaxLevelOfObjectRefinement = {{2,2,2}};

        // Mat<uint> NumSets = {{0},{0},{0}}; //First row are number of node sets, second row for number of side sets and third row for number of block sets

    };

    // Argument list for the basis functions of the B-Splines
    struct HMR_Basis_Data
    {
        BoostBitset BSplineActive;
        BoostBitset BSplineRefined;

        BoostBitset DesignBSplineActive;

        // needed for writing symmetric SDF
        BoostBitset DesignBSplineActiveSymmetric;

        // FIXME: this one is not functional and not used anywhere
        BoostBitset DesignBSplineRefined;

        BoostBitset DesignBSplineActiveLastStep;
        BoostBitset BSplineActiveLastStep;

        BoostBitset LagrangeActive;

        Mat<uint> LagrangeToBSplineMap;

        map< uint, uint > BSplineMap;

        // FIXME: these two are obsolete. Remove when fixing read_data_from_exo_file and
        //        removing read_data_from_sdf_file
        Mat<uint> NodalLocaltoGlobal;
        Mat<uint> CoordinateList;

        BoostBitset DesignLagrangeBasisActive;

        BoostBitset BSplineActiveDummy;  // temporary variable needed for MTK output
        Mat<uint> DesignBSplineActiveList;
        Mat<uint> DesignBSplineActiveListLastStep;
        Mat<uint> BSplineActiveListLastStep;
        Mat<uint> BSplineActiveList;

        map< uint, uint > DesignBSplineListMap;

        Mat<uint> LagrangeNodeProcOwner;

        // entries of T-Matrix and node IDs for each field.
        // counter starts at 1, field 0 contains length
        Mat<real> TMatrixNodalField;
        Mat<uint> IdFieldNodalField;
        Mat<real> TMatrixDesignNodalField;
        Mat<uint> IdFieldDesignNodalField;

        //Mat<uint> IdFieldDesignNodalList;

        // FIXME: this is redundant to BSplineActiveList. Figure out better way to do it
        Mat<uint> BsplineBasisList;
        uint NumberBasis;

        uint      NumberOfActiveDesignLagrangeBasis;
        uint      NumberOfNonHangingNodes = 0;
        uint      NumberOfHangingNodes = 0;

        // undocumented items below
    };

    // Argument list for for the elements of the B-Splines
    struct HMR_Mesh_Data
    {
        uint  ModelDim = 2;
        uint  Polynomial = 1;
        uint  PolynomialDesign = 1;
        uint  MaxPolynomial = 1;
        uint  NumberElements;

        Mat< uint > NumberOfBasisPerDirection;
        Mat< uint > NumberOfElementsPerDirection;

        Mat< real > Dimensions;
        Mat< real > DimensionsOffset;

        Mat<uint> NodeSet ={{0},{0}};
        Mat<uint> SideSet = {{0},{0},{0},{0},{0},{0}};
        Cell<Mat<uint>> SideSetMTK;

        // First row are number of node sets,
        // second row for number of side sets
        // and third row for number of block sets
        Mat<uint> NumSets = {{0},{0},{0}};

        uint Level = 0;
        uint LevelLastStep = 0;
        uint LevelDesignLastStep = 0;

        Mat<real> ElementLength = {{0},{0},{0}};

        uint BasisSymmetryIndex    = 0;
        uint ElementSymmetryIndex  = 0;
        uint NumberOfSymmetryMasterBasis = 0;

        Mat< real > DeleteDomain = {{0.0,0.0,0.0},{0.0,0.0,0.0}};

        // undocumented items below

        // FIXME this should be tied to either polynomial
        uint  PolynomialLagrange = 1;

        // FIXME: remove underscore and move into settings (make changes in create_mesh before)
        Mat< real > Dimensions_Orig;

    };

    // -----------------------------------------------------------------------------

    struct HMR_Proc_Data
    {
        Mat<uint> ProcCoord;
        Mat<uint> ProcNeighbor;

        // Deprecated
        //Mat<uint> ProcNeighborSendRecv;

        Mat<uint> DecompElement;
        Mat<uint> DecompBSpline;
        Mat<uint> DecompLagrangeNode;
        Mat<uint> DecompElementAura;
        Mat<uint> DecompEdgeFace;

        Mat<uint> ElementListOnProc;
        Mat<uint> ElementRefinedListOnProc;
        Mat<uint> ElementListOnProcAura;
        Mat<uint> ElementListOnProcInit;
        Mat<uint> ElementListOnProcWithAura;

        Mat<uint> ElementListOnProcLastStepDesign;
        Mat<uint> ElementListOnProcLastStepFEM;

        Mat<uint> LagrangeListOnProc;
    };

// -----------------------------------------------------------------------------

    struct HMR_Field_Data
    {
        map<uint, real> ElementDensity;
        map<uint, real> ElementVolumeDensitySDF;
        map<uint, real> ElementSurfaceDensitySDF;

        BoostBitset     NodalFieldExists;

        map<uint, real> NodalADV;
        map<uint, real> NodalPDV;
        map<uint, real> NodalLevelSet;

        Cell< map< uint, real > > NodalSDFs;

        map< uint, real > NodalSignOfSDF;

        Cell<BoostBitset> ObjectFlagSDF;
        Cell<Mat<real>>   ObjectSDF;

        // boost bitset to flag if a refinement is to be performed against this SDF
        BoostBitset       RefineAgainstObjectSDF;

        Mat<real> LagrangeADVField;
        //Mat<real> LagrangeLevelSet;

        Mat<real> BSplineCoeff;

        uint SparseSize = 0;
    };

    // -----------------------------------------------------------------------------

    class Hierarchical_Mesh_Main
    {
    protected:

    public:
        HMR_Settings         mSettings;
        HMR_Element_Data     mElementData;
        HMR_Basis_Data       mBasisData;
        HMR_Mesh_Data        mMeshData;
        HMR_Proc_Data        mProcData;
        HMR_Field_Data       mFieldData;

        //Create Object of T-Matrix
        Hierarchical_Mesh_TMatrix    mTMatrix;

        //Create Object of Refinement
        Hierarchical_Mesh_Refinement mRefinement;

        //Create Object of Face
        Hierarchical_Mesh_Face       mFace;

        //Create Object of Edge
        Hierarchical_Mesh_Edge       mEdge;

        //Create Object of Basis
        Hierarchical_Mesh_Basis      mBasis;

        //Create Object of Element
        Hierarchical_Mesh_Element    mHMRElement;

        //Create Object of sets (node, side and block set)
        Hierarchical_Mesh_MPI        mMPI;

        //Create Object of MTK
        Hierarchical_Mesh_MTK        mMTK;

        //Create Object of Output
        Hierarchical_Mesh_Output     mOutput;

        //Create Object of Input
        Hierarchical_Mesh_Input      mInput;

        //Create Object of Filter
        Hierarchical_Mesh_Filter     mFilter;

        //Create Object of Sets
        Hierarchical_Mesh_Sets       mSets;

        //Create Object of Lagrange Basis
        Lagrange_Basis               mLagrangeBasis;

        //Create Object of LagrangeElement
        Lagrange_Element             mLagrangeElement;

        //Create Object of LagrangeFace
        Lagrange_Face                mLagrangeFace;

        //Create Object of Lagrange edge
        Lagrange_Edge                mLagrangeEdge;

        //Create Object of Lagrange Filter
        Lagrange_Filter              mLagrangeFilter;

        //Create Object of BaseElement
        Base_Mesh_Element            mBaseElement;

        //Create Object of Base Edge
        Base_Mesh_Edge               mBaseEdge;

        //Create Object of Base Face
        Base_Mesh_Face               mBaseFace;

        /**
         * Hierarchical_Mesh_Main constructor
         */
        Hierarchical_Mesh_Main(
                uint const & aModelDim,
                uint const & aPolynomial,
                Mat<uint> const & aNumElements,
                Mat<real> const & aDimensions)
        {
            //Dimensions of the model
            mMeshData.ModelDim =aModelDim;

            //Polynomial degree of Bspline
            mMeshData.Polynomial = aPolynomial;

            //Polynomial degree of Lagrange basis
            mMeshData.PolynomialLagrange = aPolynomial;

            // christian: currently, design polynomial must be the same
            mMeshData.PolynomialDesign = aPolynomial;

            //Number of elements in each direction
            mMeshData.NumberOfElementsPerDirection = aNumElements;

            //Domain range of the whole problem
            mMeshData.Dimensions = aDimensions;
        };

        /**
         * Hierarchical_Mesh_Main constructor
         */
        Hierarchical_Mesh_Main()
        {
        };

        /**
         * Hierarchical_Mesh_Main destructor.
         */
        ~Hierarchical_Mesh_Main() = default;
        /**
         * Creates the Tensor mesh, based on the polynomial degree, number of elements in each direction and the domain range
         */
        void
        create_mesh();

        /**
         * Creates from a given Dimension of processors a coordinate and the neighbours
         *
         * @param[in] aDimsCartProc      Provide a vector with 3 values (Number of processors in x,y and z-direction). 0 means, that each processor returns a 0 in this coordinate direction
         *
         */
        void
        proc_coordinates_neighbors(
                Mat<uint> & aDimsCartProc);

        /**
         * Subdivides the mesh for parallel use. For each proc, lists with elements, B-Spline basis and Lagrange basis are created.
         */
        void
        mesh_decomposition();

        // Functions for base element information
        /**
         * Provides the local element ID with the local position i,j,k
         *
         * @param[in] aLevel              Level of the element.
         * @param[in] aIJKPosition        Position I,J,K.
         *
         * @param[out] Bspline            Element id.
         *
         */
        /*uint
        give_local_element_of_position(uint const & aLevel, Mat<uint> const & aIJKPosition) const
        {
            return mBaseElement.give_element_of_position(aLevel,mMeshData.ModelDim,mMeshData.NumberOfElementsPerDirectionOnProc,aIJKPosition);
        }*/

        /**
         * Provides the element with the position i,j,k
         *
         * @param[in] aLevel              Level of the element.
         * @param[in] aIJKPosition        Position I,J,K.
         *
         * @param[out] Bspline            Element id.
         *
         */
        uint
        give_element_of_position(uint const & aLevel, Mat<uint> const & aIJKPosition) const
        {
            return mBaseElement.give_element_of_position(
                    aLevel,
                    mMeshData.ModelDim,
                    mMeshData.NumberOfElementsPerDirection,
                    aIJKPosition);
        }

        /**
         * Provides the position of an element from a tensorial grid
         *
         * @param[in] aElementId            Element Id.
         *
         * @param[out] Element_position        Position I,J,K.
         *
         */
        Mat<uint>
        give_position_of_element(uint const & aElementId) const
        {
            return mBaseElement.give_position_of_element(
                    aElementId,
                    mMeshData.ModelDim,
                    mMeshData.NumberOfElementsPerDirection);
        }

        /**
         * Provides the level of the element from the hierarchical mesh.
         *
         * @param[in] aElementId            Element Id.
         *
         * @param[out] level        Level of the elemen Id t from the hierarchical mesh.
         *
         */
        uint
        give_element_level(uint const & aElementId) const
        {
            return mBaseElement.give_element_level(
                    aElementId,
                    mMeshData.ModelDim,
                    mMeshData.NumberOfElementsPerDirection);
        }

        /**
         * Provides the parent of an element
         *
         * @param[in] aElementId            Element Id.
         *
         * @param[out] Parent        Parent of an element Id.
         *
         */
        uint
        give_parent_of_element(uint const & aElementId)
        {
            return  mBaseElement.give_parent_of_element(
                    aElementId,
                    mMeshData.ModelDim,
                    mMeshData.NumberOfElementsPerDirection);
        }

        /**
         * Provides the parent of an element for a specific level
         *
         * @param[in] aElementId            Element Id.
         * @param[in] aWhichLevel         User defined level on which the parent needs to be determined
         *
         * @param[out] Parent        Parent of an element.
         *
         */
        uint
        give_parent_of_level_x(uint const & aElementId, uint const & aWhichLevel)
        {
            return mBaseElement.give_parent_of_level_x(
                    aElementId,
                    mMeshData.ModelDim,
                    mMeshData.NumberOfElementsPerDirection,
                    aWhichLevel);
        }

        /**
         * Provides the parent of an element and the relation to the child
         *
         * @param[in] aElementId            Element Id.
         *
         * @param[out] Parent          A vector with two entries. First entry is the parent of the element Id and second entry is the child of the parent relation (Which child is element Id of parent)
         *
         */
        Mat<uint>
        give_parent_child_realation_of_element(uint const & aElementId)
        {
            return  mBaseElement.give_parent_child_realation_of_element(
                    aElementId,
                    mMeshData.ModelDim,
                    mMeshData.NumberOfElementsPerDirection);
        }

        /**
         * Provides the number of elements within all levels from level 0 until "aLevel"
         *
         * @param[in] aLevel              Level of the element.
         *
         * @param[out] element_number        Number of elements within all Levels until aLevel.
         *
         */
        uint
        give_number_of_elements(uint const & aLevel) const
        {
            return mBaseElement.give_number_of_elements(
                    aLevel,
                    mMeshData.ModelDim,
                    mMeshData.NumberOfElementsPerDirection);
        }

        /**
         * Provides the children of an element
         *
         * @param[in] aElementId            Element Id.
         *
         * @param[out] Children        Childrens of an element (4 children for Dim=2 and 8 children for Dim=3).
         *
         */
        Mat<uint>
        give_children_of_element(uint const & aElementId) const
        {
            return mBaseElement.give_children_of_element(
                    aElementId,
                    mMeshData.ModelDim,
                    mMeshData.NumberOfElementsPerDirection);
        }

        /**
         * Provides the neighbors of an element
         *
         * @param[in] aElementId            Element Id.
         * @param[in] aBuffer             Provides the number of layers around the element
         *
         * @param[out] Element_neighbour        Element plus neighbors (9 for aModelDim = 2, 27 for aModelDim = 3) if a buffer of one is used.\n
         *                                2D:...----------------------------------...3D:....-----plain 1----------------------.....-----plain 2-------------------------\n.....-----plain 3-------------------------\n
         *                                ......|Element(6) Element(7) Element(8)|..........|Element(6) Element(7) Element(8)|.....|Element(15) Element(16) Element(17)|\n.....|Element(24) Element(25) Element(26)|\n
         *                                ......|Element(3) Element(4) Element(5)|..........|Element(3) Element(4) Element(5)|.....|Element(12) Element(13) Element(14)|\n.....|Element(21) Element(22) Element(23)|\n
         *                                ......|Element(0) Element(1) Element(2)|..........|Element(0) Element(1) Element(2)|.....|Element(9) Element(10) Element(11) |\n.....|Element(18) Element(19) Element(20)|\n
         *                                ......|--------------------------------|..........|--------------------------------|.....|-----------------------------------|  .....|-----------------------------------|\n
         *                                2D: Element(4) = aElementId, 3D: Element(13) = aElementId;
         *
         */
        Mat<uint>
        give_neighbor_of_element(const uint & aElementId, const uint & aBuffer)
        {
            return mBaseElement.give_neighbor_of_element(
                    aElementId,
                    mMeshData.ModelDim,
                    aBuffer,
                    mMeshData.NumberOfElementsPerDirection);
        }

        /**
         * Provides the ownership of an element
         *
         * @param[in] aElementId            Element Id.
         *
         * @param[out] Rank        Proc rank
         *
         */
        /*uint
        give_element_owner(uint const & aElementId)
        {
            return mBaseElement.give_element_owner(
                    aElementId,
                    mMeshData.ModelDim,
                    mMeshData.NumberOfElementsPerDirection,
                    mProcData.DecompElement,
                    mProcData.ProcNeighborSendRecv);
        } */

        /**
         * Provides a list, which share an element (owner plus aura procs) - An proc can identify owner and share only if the element lives on the proc or in the aura of that proc
         *
         * @param[in] aElementId    Element Id.
         *
         * @param[out] Share        Vector with procs, which share the element
         *
         */
        /* Mat<uint>
        give_element_share(uint const & aElementId)
        {
            return mBaseElement.give_element_share(
                    aElementId,
                    mMeshData.ModelDim,
                    mMeshData.NumberOfElementsPerDirection,
                    mProcData.DecompElement,
                    mProcData.ProcNeighborSendRecv);
        } */

        /**
         * Provides the neighbours of an element, which are connected through the faces. This function uses the function "give_neighbour_of_element" and extracts only the nessecary elements
         *
         * @param[in] aElementId            Element number.
         *
         * @param[out] Element_neighbour        Element plus neighbours (4 for aModelDim = 2, 6 for aModelDim = 3).
         *
         */
        Mat<uint>
        give_face_neighbor_of_element(uint const & aElementId) const
        {
            return mBaseElement.give_face_neighbor(
                    aElementId,
                    mMeshData.ModelDim,
                    mMeshData.NumberOfElementsPerDirection);
        }

        //-----------------------------------------------------------------------------------
        //Element functionalities with respect to HMR
        /**
         * Provides the basis function IDs of an element from a tensor grid
         *
         * @param[in] aElementId            Element Id.
         *
         * @param[out] basis        Basis function IDs of an element.
         *
         */
        Mat<uint>
        give_basis_of_element(uint const & aElementId) const
        {
            return mHMRElement.give_basis_of_element(
                    aElementId,
                    mMeshData.ModelDim,
                    mMeshData.Polynomial,
                    mMeshData.NumberOfElementsPerDirection);
        }

        /**
         * Provides all active face neighbors (includes all levels)
         *
         * @param[in] aElementId            Element number.
         *
         * @param[out] Element_neighbor   Active neighbors
         *
         */
        Mat<uint>
        give_active_face_neighbor_of_element(uint const & aElementId) const
        {
            return mHMRElement.give_active_face_neighbor_of_element(
                    aElementId,
                    mMeshData.ModelDim,
                    mMeshData.NumberOfElementsPerDirection,
                    mMeshData.Level,
                    mElementData.ElementActive);
        }

        /**
         * Provides all active neighbors (includes all levels) (includes also neighbors conntected through nodes
         *
         * @param[in] aElementId            Element number.
         *
         * @param[out] Element_neighbor   Active neighbors
         *
         */
        Mat<uint>
        give_active_neighbor_of_element(uint const & aElementId) const
        {
            return mHMRElement.give_active_neighbor_of_element(
                    aElementId,
                    mMeshData.ModelDim,
                    mMeshData.NumberOfElementsPerDirection,
                    mMeshData.Level,
                    mElementData.ElementActive);
        }

        /**
         * Provides the coordinate in the middle of an element, it is needed to have a rough idea where the element sits, needed for "Delete Domain in Hierarchical_Mesh_Main::create_mesh()"
         *
         * @param[in] aElementId        A specific element.
         *
         * @param[out] coordinates       Coordinates in x,y,z direction.
         *
         */
        Mat<real>
        give_middlecoordinate_from_element(uint const & aElementId)
        {
            return mHMRElement.give_middlecoordinate_from_element(
                    aElementId,
                    mMeshData.ModelDim,
                    mMeshData.NumberOfElementsPerDirection,
                    mMeshData.Dimensions,
                    mMeshData.DimensionsOffset);
        }

        /**
         * Provides a list of active elements
         *
         * @param[out] ElementListOnProc      Updated list of elements, which are active on the current proc
         *
         */
        Mat<uint>
        give_active_elements() const
        {
            return mRefinement.give_active_elements(
                    mMeshData.ModelDim,
                    mMeshData.Polynomial,
                    mMeshData.NumberOfElementsPerDirection,
                    mMeshData.Level,
                    mProcData.ElementListOnProcInit,
                    mElementData.ElementActive);
        }

        /**
         * Provides a list of active elements with the aura
         *
         * @param[out] ElementListOnProc      Updated list of elements, which are active on the current proc with the aura
         *
         */
        Mat<uint>
        give_active_elements_with_aura() const
        {
            return mRefinement.give_active_elements(
                    mMeshData.ModelDim,
                    mMeshData.Polynomial,
                    mMeshData.NumberOfElementsPerDirection,
                    mMeshData.Level,
                    mProcData.ElementListOnProcWithAura,
                    mElementData.ElementActive);
        }

        /**
         * Provides a list of passive elements (uses the same function as for active elements, but the passive bitset is used!)
         *
         *
         * @param[out] ElementListOnProc      Updated list of elements, which are active on the current proc
         *
         * @TODO change name of this function to give_refined_elements
         */
        Mat<uint>
        give_passive_elements() const
        {
            return mRefinement.give_active_elements(
                    mMeshData.ModelDim,
                    mMeshData.Polynomial,
                    mMeshData.NumberOfElementsPerDirection,
                    mMeshData.Level,
                    mProcData.ElementListOnProcInit,
                    mElementData.ElementRefined);
        }
        // Functions for basis information-----------------------------------------------------------------
        /**
         * Provides the number of basis functions within all levels from level 0 until aLevel
         *
         * @param[in] aLevel              Level of the basis functions.
         *
         * @param[out] basis_number        Number of basis functions within all Levels until aLevel.
         *
         */
        uint
        give_number_of_basis(const uint & aLevel)
        {
            return mBasis.give_number_of_basis(
                    aLevel,
                    mMeshData.ModelDim,
                    mMeshData.Polynomial,                        // FIXME: must be maxPolynomial ?
                    mMeshData.NumberOfElementsPerDirection);
        }

        /**
         * Provides the neighbours of an basis
         *
         * @param[in] aBasisId              Basis function Id.
         * @param[in] aBuffer             Provides the number of layers around the basis
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
        give_neighbor_of_basis(const uint & aBasisId, const uint & aBuffer)
        {
            return mBasis.give_neighbor_of_basis(
                    aBasisId,
                    mMeshData.ModelDim,
                    mMeshData.Polynomial,
                    aBuffer,
                    mMeshData.NumberOfElementsPerDirection);
        }

        /**
         * Provides the level of the basis function
         *
         * @param[in] aBasisId            Basis function Id.
         *
         * @param[out] level        Level of the element from the hierarchical mesh.
         *
         */
        uint
        give_basis_level(const uint & aBasisId) const
        {
            return mBasis.give_basis_level(
                    aBasisId,
                    mMeshData.ModelDim,
                    mMeshData.Polynomial,
                    mMeshData.NumberOfElementsPerDirection);
        }

        /**
         * Provides the position of a basis function from a tensorial grid
         *
         * @param[in] aBasisId            Basis function Id.
         *
         * @param[out] Basis_position        Position I,J,K.
         *
         */
        Mat<uint>
        give_position_of_basis(const uint & aBasisId) const
        {
            return mBasis.give_position_of_basis(
                    aBasisId,
                    mMeshData.ModelDim,
                    mMeshData.Polynomial,
                    mMeshData.NumberOfElementsPerDirection);
        }

        /**
         * Provides the basis function Id with the position i,j,k
         *
         * @param[in] aLevel             Level of the basis function Id.
         * @param[in] aIJKPosition        Position I,J,K.
         *
         * @param[out] tBasis       Basis function Id.
         *
         */
        uint
        give_basis_of_position( const uint & aLevel, const Mat<uint> & aIJKPosition) const
        {
            return mBasis.give_basis_of_position(
                    aLevel,
                    mMeshData.ModelDim,
                    mMeshData.Polynomial,  // FIXME Max Polynomial ?
                    mMeshData.NumberOfElementsPerDirection,
                    aIJKPosition);
        }

        /**
         * Provides the basis function Id of the parent (Works only for a linear polynomial degree) (Returns UINT_MAX if there is no parent)
         *
         * @param[in] aBasisId            Basis function number.
         *
         * @param[out] tBasis            Basis function Id of the parent.
         *
         */
        uint
        give_basis_of_parent( const uint & aBasisId )
        {
            return mBasis.give_basis_of_parent(
                    aBasisId,
                    mMeshData.ModelDim,
                    mMeshData.Polynomial,
                    mMeshData.NumberOfElementsPerDirection);
        }

        /**
         * Provides the elements, which have support with the basis function from a tensorial grid
         *
         * @param[in] aBasisId            Basis function Id.
         *
         * @param[out] tElements        Elements with support of basis function Id.
         *
         */
        Mat<uint>
        give_element_of_basis( const uint & aBasisId ) const
        {
            return mBasis.give_element_of_basis(
                    aBasisId,
                    mMeshData.ModelDim,
                    mMeshData.Polynomial,
                    mMeshData.NumberOfElementsPerDirection);
        }

        /**
         * Provides the coordinates for a specific basis function Id
         *
         * @param[in] aBasisId                Basis function Id
         *
         * @param[out] tCoordinates       Coordinates in x,y,z direction.
         *
         */
        Mat<real>
        give_coordinate_from_basis( const uint & aBasisId ) const
        {
            return mBasis.give_coordinate_from_basis(
                    aBasisId,
                    mMeshData.ModelDim,
                    mMeshData.Polynomial,
                    mMeshData.NumberOfElementsPerDirection,
                    mMeshData.Dimensions,
                    mMeshData.DimensionsOffset);
        }

        /**
         * Defines owner of a basis function Id
         *
         * @param[in] aBasisId                  Basis function Id
         *
         * @param[out] tProcOwner           Owner of the basis function id
         *
         */
        uint
        give_basis_owner( const uint & aBasisId) const
        {
            return mBasis.give_basis_proc_owner(
                    aBasisId,
                    mMeshData.ModelDim,
                    mMeshData.Polynomial,
                    mMeshData.NumberOfElementsPerDirection,
                    mProcData.ProcNeighbor,
                    mProcData.DecompBSpline);
        }

        /**
         * Defines a list of processors, which share a basis function Id
         *
         * @param[in] aBasisId                 Basis function Id

         *
         * @param[out] tProcShare           Who shares this basis function
         *
         */
        Mat<uint>
        give_basis_share( const  uint & aBasisId ) const
        {
            return mBasis.give_basis_share(
                    aBasisId,
                    mMeshData.ModelDim,
                    mMeshData.Polynomial,
                    mMeshData.NumberOfElementsPerDirection,
                    mProcData.ProcNeighbor,
                    mProcData.DecompBSpline);
        }

        //---------------------------------------------------------------------------
        //Lagrange element functionalities
        /**
         * Provides the shape function IDs of an element from a tensorial grid
         *
         * @param[in] aElementId            Element Id.
         *
         * @param[out] basis        Basis function IDs of an element.
         *
         */
        Mat<uint>
        give_lagrange_basis_of_element( const uint & aElementId ) const
        {
            return mLagrangeElement.give_basis_of_element(
                    aElementId,
                    mMeshData.ModelDim,
                    mMeshData.Polynomial,
                    mMeshData.NumberOfElementsPerDirection);
        }
        /**
         * Provides the "active" shape function IDs of an element from a tensorial grid (They are not really active, but always the lowest possible parent is provided )
         *
         * @param[in] aElementId            Element Id.
         *
         * @param[out] basis        Basis function IDs of an element.
         *
         */
        Mat<uint>
        give_active_lagrange_basis_of_element( const uint & aElementId ) const;

        /**
         * Provides the local "active" shape function IDs of an element from a tensorial grid (They are not really active, but always the lowest possible parent is provided )
         *
         * @param[in] aElementId            Element Id.
         *
         * @param[out] basis        Basis function IDs of an element.
         */
        //Mat<uint>
        //give_local_active_lagrange_basis_of_element( const uint & aElementId ) const;

        //---------------------------------------------------------------------------
        //Lagrange basis functionalities
        /**
         * Provides the number of lagrange basis functions within all levels from level 0 until aLevel
         *
         * @param[in] aLevel              Level of the basis functions.
         *
         * @param[out] basis_number        Number of lagrange basis functions within all Levels until aLevel.
         *
         */
        uint
        give_number_of_lagrange_basis( const uint & aLevel ) const
        {
            return mLagrangeBasis.give_number_of_basis(
                    aLevel,
                    mMeshData.ModelDim,
                    mMeshData.Polynomial,
                    mMeshData.NumberOfElementsPerDirection);
        }

        /**
         * Provides the level of the lagrange basis function
         *
         * @param[in] aBasisId      Lagrange basis function Id.
         *
         * @param[out] level        Level of the element from the hierarchical mesh.
         *
         */
        uint
        give_lagrange_basis_level(const uint & aBasisId ) const
        {
            return mLagrangeBasis.give_basis_level(
                    aBasisId,
                    mMeshData.ModelDim,
                    mMeshData.Polynomial,
                    mMeshData.NumberOfElementsPerDirection );
        }

        /**
         * Provides the local level of the lagrange basis function
         *
         * @param[in] aBasisId      Lagrange basis function Id.
         *
         * @param[out] level        Level of the element from the hierarchical mesh.
         *
         */
        /* uint
        give_local_lagrange_basis_level( const uint & aBasisId ) const
        {
            return mLagrangeBasis.give_basis_level(
                    aBasisId,
                    mMeshData.ModelDim,
                    mMeshData.Polynomial,
                    mMeshData.NumberOfElementsPerDirectionOnProc );
        } */

        /**
         * Provides the position of a lagrange basis function from a tensorial grid
         *
         * @param[in] aBasisId            Lagrange basis function Id.
         *
         * @param[out] Basis_position        Position I,J,K.
         *
         */
        Mat<uint>
        give_position_of_lagrange_basis( const  uint & aBasisId ) const
        {
            return mLagrangeBasis.give_position_of_basis(
                    aBasisId,
                    mMeshData.ModelDim,
                    mMeshData.Polynomial,
                    mMeshData.NumberOfElementsPerDirection );
        }
        /**
         * Provides the lagrange basis function Id with the position i,j,k
         *
         * @param[in] aLevel             Level of the lagrange basis function Id.
         * @param[in] aIJKPosition        Position I,J,K.
         *
         * @param[out] tBasis       Lagrange basis function Id.
         *
         */
        uint
        give_lagrange_basis_of_position( const uint & aLevel, const Mat<uint> & aIJKPosition )
        {
            return mLagrangeBasis.give_basis_of_position(
                    aLevel,
                    mMeshData.ModelDim,
                    mMeshData.Polynomial,
                    mMeshData.NumberOfElementsPerDirection,
                    aIJKPosition);
        }

        /**
         * Provides the local lagrange basis function Id with the position i,j,k
         *
         * @param[in] aLevel             Level of the lagrange basis function Id.
         * @param[in] aIJKPosition        Position I,J,K.
         *
         * @param[out] tBasis       Lagrange basis function Id.
         *
         */
        /* uint
        give_local_lagrange_basis_of_position( const uint & aLevel, const Mat<uint> & aIJKPosition ) const
        {
            return mLagrangeBasis.give_basis_of_position(
                    aLevel,
                    mMeshData.ModelDim,
                    mMeshData.Polynomial,
                    mMeshData.NumberOfElementsPerDirectionOnProc,
                    aIJKPosition);
        } */

        /**
         * Provides the lagrange basis function Id of the parent (Works only for a linear polynomial degree) (Returns UINT_MAX if there is no parent)
         *
         * @param[in] aBasisId            Basis function number.
         *
         * @param[out] tBasis            Lagrange basis function Id of the parent.
         *
         */
        uint
        give_lagrange_basis_of_parent( const uint & aBasisId ) const
        {
            return mLagrangeBasis.give_basis_of_parent(
                    aBasisId,
                    mMeshData.ModelDim,
                    mMeshData.Polynomial,
                    mMeshData.NumberOfElementsPerDirection);
        }

        /**
         * Provides the local lagrange basis function Id of the parent (Works only for a linear polynomial degree) (Returns UINT_MAX if there is no parent)
         *
         * @param[in] aBasisId            Basis function number.
         *
         * @param[out] tBasis            Lagrange basis function Id of the parent.
         *
         */
        /* uint
        give_local_lagrange_basis_of_parent( const uint & aBasisId ) const
        {
            return mLagrangeBasis.give_basis_of_parent(
                    aBasisId,
                    mMeshData.ModelDim,
                    mMeshData.Polynomial,
                    mMeshData.NumberOfElementsPerDirectionOnProc);
        } */
        /**
         * Provides for a list of lagrange basis function Id the parent basis ids (if possible)
         *
         * @param[in] aBasisList         List of basis functions
         *
         * @param[out] tBasisList        List of updated basis functions
         *
         */
        Mat<uint>
        give_parents_of_lagrange_basis( const Mat<uint> & aBasisList ) const;

        /**
         * Provides for a list of lagrange basis function Id the parent basis ids (if possible)
         *
         * @param[in] aBasisList         List of basis functions
         *
         * @param[out] tBasisList        List of updated basis functions
         *
         */
        // Mat<uint>
        // give_local_parents_of_lagrange_basis( const Mat<uint> & aBasisList ) const;

        /**
         * Provides the neighbours of an basis
         *
         * @param[in] aBasisId              Basis function Id.
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
        give_neighbor_of_lagrange_basis( const uint & aBasisId, const uint & aBuffer )
        {
            return mLagrangeBasis.give_neighbor_of_basis(
                    aBasisId,
                    mMeshData.ModelDim,
                    mMeshData.Polynomial,
                    aBuffer,
                    mMeshData.NumberOfElementsPerDirection);
        }

        /**
         * Provides the elements, which have support with the basis function from a tensorial grid
         *
         * @param[in] aBasisId            Basis function Id.
         *
         * @param[out] tElements        Elements with support of basis function Id.
         *
         */
        Mat<uint>
        give_element_of_lagrange_basis( const uint & aBasisId ) const
        {
            return mLagrangeBasis.give_element_of_basis(
                    aBasisId,
                    mMeshData.ModelDim,
                    mMeshData.Polynomial,
                    mMeshData.NumberOfElementsPerDirection);
        }

        /**
         * Provides the coordinates for a specific basis function Id
         *
         * @param[in] aBasisId                Basis function Id
         *
         * @param[out] tCoordinates       Coordinates in x,y,z direction.
         *
         */
        Mat<real>
        give_coordinate_from_lagrange_basis( const uint & aBasisId ) const
        {
            return mLagrangeBasis.give_coordinate_from_basis(aBasisId,
                    mMeshData.ModelDim,
                    mMeshData.Polynomial,
                    mMeshData.NumberOfElementsPerDirection,
                    mMeshData.Dimensions,
                    mMeshData.DimensionsOffset);
        }

        /**
         * Defines owner of a lagrange basis function Id
         *
         * @param[in] aBasisId                  Lagrange basis function Id
         *
         * @param[out] tProcOwner           Owner of the basis function id
         *
         */
        uint
        give_lagrange_basis_owner( const uint & aBasisId ) const
        {
            return mLagrangeBasis.give_basis_proc_owner(
                    aBasisId,
                    mMeshData.ModelDim,
                    mMeshData.Polynomial,
                    mMeshData.NumberOfElementsPerDirection,
                    mProcData.ProcNeighbor,
                    mProcData.DecompLagrangeNode);
        }

        /**
         * Defines a list of processors, which share a lagrange basis function Id
         *
         * @param[in] aBasisId                 Lagrange basis function Id

         *
         * @param[out] tProcShare           Who shares this basis function
         *
         */
        Mat<uint>
        give_lagrange_basis_share( const uint & aBasisId ) const
        {
            return mLagrangeBasis.give_basis_share(
                    aBasisId,
                    mMeshData.ModelDim,
                    mMeshData.Polynomial,
                    mMeshData.NumberOfElementsPerDirection,
                    mProcData.ProcNeighbor,
                    mProcData.DecompLagrangeNode);
        }

        // Functions for edge information-----------------------------------------------------------------------------------------------------------
        /**
         * Provides the number of edges in x-direction
         *
         * @param[in] aLevel              Level of the element.
         *
         * @param[out] tEdgeNumber        Number of edges in x-direction.
         *
         */
        uint
        give_number_of_edges_x( const uint & aLevel ) const
        {
            return mBaseEdge.give_number_of_edges_x(
                    aLevel,
                    mMeshData.ModelDim,
                    mMeshData.NumberOfElementsPerDirection);
        }

        /**
         * Provides the number of edges in y-direction
         *
         * @param[in] aLevel              Level of the element.
         *
         * @param[out] tEdgeNumber        Number of edges in y-direction.
         *
         */
        uint
        give_number_of_edges_y( const uint & aLevel ) const
        {
            return mBaseEdge.give_number_of_edges_y(
                    aLevel,
                    mMeshData.ModelDim,
                    mMeshData.NumberOfElementsPerDirection);
        }

        /**
         *  Provides the number of edges in z-direction
         *
         * @param[in] aLevel              Level of the element.
         *
         * @param[out] tEdgeNumber        Number of edges in z-direction.
         *
         */
        uint
        give_number_of_edges_z( const uint & aLevel ) const
        {
            return mBaseEdge.give_number_of_edges_z(
                    aLevel,
                    mMeshData.ModelDim,
                    mMeshData.NumberOfElementsPerDirection);
        }

        /**
         * Provides the number of edges within all levels from level 0 until aLevel
         *
         * @param[in] aLevel              Level of the basis functions.
         *
         * @param[out] edges_number        Number of edges within all Levels until aLevel.
         *
         */
        uint
        give_number_of_edges( const uint & aLevel) const
        {
            return mBaseEdge.give_number_of_edges(
                    aLevel,
                    mMeshData.ModelDim,
                    mMeshData.NumberOfElementsPerDirection);
        }

        /**
         * Provides the edge number in x-direction of the position i,j,k
         *
         * @param[in] aLevel              Level of the element.
         * @param[in] aIJKPosition        Position I,J,K.
         *
         * @param[out] tEdgeNumber        Number of edge in x-direction.
         *
         */
        uint
        give_edge_x_of_position( const uint & aLevel, const Mat<uint> & aIJKPosition ) const
        {
            return mBaseEdge.give_edge_x_of_position(
                    aLevel,
                    mMeshData.ModelDim,
                    mMeshData.NumberOfElementsPerDirection,
                    aIJKPosition);
        }

        /**
         * Provides the edge number in y-direction of the position i,j,k
         *
         * @param[in] aLevel              Level of the element.
         * @param[in] aIJKPosition        Position I,J,K.
         *
         * @param[out] tEdgeNumber        Number of edge in y-direction.
         *
         */
        uint
        give_edge_y_of_position( const uint & aLevel, const Mat<uint> & aIJKPosition ) const
        {
            return mBaseEdge.give_edge_y_of_position(
                    aLevel,
                    mMeshData.ModelDim,
                    mMeshData.NumberOfElementsPerDirection,
                    aIJKPosition);
        }

        /**
         * Provides the edge number in z-direction of the position i,j,k
         *
         * @param[in] aLevel              Level of the element.
         * @param[in] aIJKPosition        Position I,J,K.
         *
         * @param[out] tEdgeNumber        Number of edge in z-direction.
         *
         */
        uint
        give_edge_z_of_position( const uint & aLevel, const Mat<uint> & aIJKPosition ) const
        {
            return mBaseEdge.give_edge_z_of_position(
                    aLevel,
                    mMeshData.ModelDim,
                    mMeshData.NumberOfElementsPerDirection,
                    aIJKPosition);
        }

        /**
         * Provides the level of the edge from the hierarchical mesh.
         *
         * @param[in] aEdgeId            Edge number.
         *
         * @param[out] level        Level of the edge from the hierarchical mesh.
         *
         */
        uint
        give_edge_level( const uint & aEdgeId ) const
        {
            return mBaseEdge.give_edge_level(
                    aEdgeId,
                    mMeshData.ModelDim,
                    mMeshData.NumberOfElementsPerDirection);
        }

        /**
         * Provides the edges, which are connected to an element.
         *
         * @param[in] aElement            Element number.
         *
         * @param[out] Edges        Edges, which are connected to an element. In 2D: Edges(4,1) = [left edge_X, right edge_X, bottom edge_Y, top_edge_Y], In 2D: Edges(12,1) = [bottom left edge_X, bottom right edge_X, top left edge_X, top right edge_X, front bottom edge_Y, front top_edge_Y, back bottom edge_Y, back top_edge_Y, frontleft edge_z, front right edge_z, back left edge_z, back right edge_z]
         *
         */
        Mat<uint>
        give_element_edges(  const uint & aElementId ) const
        {
            return mBaseEdge.give_element_edges(
                    aElementId,
                    mMeshData.ModelDim,
                    mMeshData.NumberOfElementsPerDirection);
        }

        /**
         * Provides the active elements, which are connected to an edge.
         *
         * @param[in] aEdgeId            Edge Id.
         * @param[in] aModelDim                Number of dimensions.
         * @param[in] aNumberOfElementsPerDirection        NumElements=[Number of elements in x,y,z direction].
         *
         * @param[out] Element        Elements, which are connected to an edge.
         *
         */
        Mat<uint>
        give_active_elements_of_edge( const uint & aEdgeId ) const;

        /**
         * Provides the position of an edge
         *
         * @param[in] aEdgeId            Edge Id.
         *
         * @param[out] Position         Position I,J,K.
         *
         */
        Mat<uint>
        give_edge_position( const uint & aEdgeId ) const
        {
            return mBaseEdge.give_edge_position(
                    aEdgeId,
                    mMeshData.ModelDim,
                    mMeshData.NumberOfElementsPerDirection);
        }

        /**
         * Provides the ownership of an edge
         *
         * @param[in] aEdgeId            Edge number.
         *
         * @param[out] Rank        Proc rank
         *
         */
        uint
        give_edge_owner(const uint & aEdgeId) const
        {
            return mBaseEdge.give_edge_owner(
                    aEdgeId,
                    mMeshData.ModelDim,
                    mMeshData.NumberOfElementsPerDirection,
                    mProcData.DecompEdgeFace,
                    mProcData.ProcNeighbor);
        }

        /**
         * Provides the a list, which share an edge
         *
         * @param[in] aEdgeId            Edge number.
         *
         * @param[out] Share        Vector with procs, which share the edge
         *
         */
        Mat<uint>
        give_edge_share( const uint & aEdgeId ) const
        {
            return mBaseEdge.give_edge_share(
                    aEdgeId,
                    mMeshData.ModelDim,
                    mMeshData.NumberOfElementsPerDirection,
                    mProcData.DecompEdgeFace,
                    mProcData.ProcNeighbor);
        }

        /**
         * Provides the elements, which are connected to an edge.
         *
         * @param[in] aEdgeId            Edge Id.
         * @param[in] aModelDim                Number of dimensions.
         * @param[in] aNumberOfElementsPerDirection        NumElements=[Number of elements in x,y,z direction].
         *
         * @param[out] Element        Elements, which are connected to an edge.
         *
         */
        Mat<uint>
        give_elements_of_edge(  const uint & aEdgeId ) const
        {
            return mBaseEdge.give_elements_of_edge(
                    aEdgeId,
                    mMeshData.ModelDim,
                    mMeshData.NumberOfElementsPerDirection);
        }

        //Functions for B-spline basis and edge---------------------------------
        /**
         * Provides the nodes of an edge
         *
         * @param[in] aEdgeId            Edge number.
         *
         * @param[out] Nodes         Node id's, which are connected to an edge (Only for linear possible)
         *
         */
        Mat<uint>
        give_edge_nodes( const uint & aEdgeId ) const
        {
            return mEdge.give_edge_nodes(
                    aEdgeId,
                    mMeshData.ModelDim,
                    mMeshData.NumberOfElementsPerDirection);
        }

        //Functions for Lagrange basis and edge---------------------------------
        /**
         * Provides the nodes of an edge
         *
         * @param[in] aEdgeId            Edge number.
         *
         * @param[out] Nodes         Node id's, which are connected to an edge (Only for linear possible)
         *
         */
        Mat<uint>
        give_edge_lagrange_nodes( const uint & aEdgeId ) const
        {
            return mLagrangeEdge.give_edge_nodes(
                    aEdgeId,
                    mMeshData.ModelDim,
                    mMeshData.Polynomial,
                    mMeshData.NumberOfElementsPerDirection);
        }

        /**
         * Provides the "active" nodes of an edge (They are not really active, but always the lowest possible parent is provided )
         *
         * @param[in] aEdgeId            Edge number.
         *
         * @param[out] Nodes         Node id's, which are connected to an edge (Only for linear possible)
         *
         */
        Mat<uint>
        give_edge_active_lagrange_nodes( const uint & aEdgeId) const;

        //Functions for face information----------------------------------------
        /**
         * Provides the number of faces in x-direction
         *
         * @param[in] aLevel              Level of the element.
         *
         * @param[out] Face_number        Number of faces in x-direction.
         *
         */
        uint
        give_number_of_faces_x(uint & aLevel)
        {
            return mBaseFace.give_number_of_faces_x(aLevel,mMeshData.ModelDim,mMeshData.NumberOfElementsPerDirection);
        }

        /**
         * Provides the number of faces in y-direction
         *
         * @param[in] aLevel              Level of the element.
         *
         * @param[out] Face_number        Number of faces in y-direction.
         *
         */
        uint
        give_number_of_faces_y(uint & aLevel)
        {
            return mBaseFace.give_number_of_faces_y(aLevel,mMeshData.ModelDim,mMeshData.NumberOfElementsPerDirection);
        }

        /**
         * Provides the number of faces in z-direction
         *
         * @param[in] aLevel              Level of the element.
         *
         * @param[out] Face_number        Number of faces in z-direction.
         *
         */
        uint
        give_number_of_faces_z(uint & aLevel)
        {
            return mBaseFace.give_number_of_faces_z(aLevel,mMeshData.ModelDim,mMeshData.NumberOfElementsPerDirection);
        }

        /**
         * Provides the number of faces
         *
         * @param[in] aLevel              Level of the element.
         *
         * @param[out] Face_number        Number of faces in z-direction.
         *
         */
        uint
        give_number_of_faces(uint & aLevel)
        {
            return mBaseFace.give_number_of_faces(aLevel,mMeshData.ModelDim,mMeshData.NumberOfElementsPerDirection);
        }

        /**
         * Provides the level of the face from the hierarchical mesh.
         *
         * @param[in] aFaceId            Face number.
         *
         * @param[out] level        Level of the face from the hierarchical mesh.
         *
         */
        uint
        give_face_level(uint & aEdge)
        {
            return mBaseFace.give_face_level(aEdge,mMeshData.ModelDim,mMeshData.NumberOfElementsPerDirection);
        }

        /**
         * Provides the face number in x-direction of the position i,j,k
         *
         * @param[in] aLevel              Level of the element.
         * @param[in] aIJKPosition        Position I,J,K.
         *
         * @param[out] Face_number        Number of face in x-direction.
         *
         */
        uint
        give_face_x_of_position(uint const & aLevel, Mat<uint> const & aIJKPosition)
        {
            return mBaseFace.give_face_x_of_position(aLevel,mMeshData.ModelDim,mMeshData.NumberOfElementsPerDirection,aIJKPosition);
        }

        /**
         * Provides the face number in y-direction of the position i,j,k
         *
         * @param[in] aLevel              Level of the element.
         * @param[in] aIJKPosition        Position I,J,K.
         *
         * @param[out] Face_number        Number of face in y-direction.
         *
         */
        uint
        give_face_y_of_position(uint const & aLevel, Mat<uint> const & aIJKPosition)
        {
            return mBaseFace.give_face_y_of_position(aLevel,mMeshData.ModelDim,mMeshData.NumberOfElementsPerDirection,aIJKPosition);
        }

        /**
         * Provides the face number in z-direction of the position i,j,k
         *
         * @param[in] aLevel              Level of the element.
         * @param[in] aIJKPosition        Position I,J,K.
         *
         * @param[out] Face_number        Number of face in z-direction.
         *
         */
        uint
        give_face_z_of_position(uint const & aLevel, Mat<uint> const & aIJKPosition)
        {
            return mBaseFace.give_face_z_of_position(aLevel,mMeshData.ModelDim,mMeshData.NumberOfElementsPerDirection,aIJKPosition);
        }

        /**
         * Provides the faces, which are connected to an element.
         *
         * @param[in] aElement            Element number.
         *
         * @param[out] Faces        Faces, which are connected to an element. In 2D: Faces(4,1) = [left face_X, right face_X, bottom face_Y, top_face_Y], In 2D: faces(6,1) = [left face_X, right face_X,  bottom face_Y, top_face_Y, front face_z, back face_z]
         *
         */
        Mat<uint>
        give_element_faces(uint const & aElementId) const
        {
            return mBaseFace.give_element_faces(aElementId,mMeshData.ModelDim,mMeshData.NumberOfElementsPerDirection);
        }

        /**
         * Provides the elements, which are connected to an edge.
         *
         * @param[in] aEdge            Edge Id.
         *
         * @param[out] Element        Elements, which are connected to an edge.
         *
         */
        Mat<uint>
        give_elements_of_face(uint const & aFaceId) const
        {
            return mBaseFace.give_elements_of_face(aFaceId,mMeshData.ModelDim,mMeshData.NumberOfElementsPerDirection);
        }

        /**
         * Provides the active elements, which are connected to an edge.
         *
         * @param[in] aEdge            Edge Id.
         *
         * @param[out] Element        Elements, which are connected to an edge.
         *
         */
        Mat<uint>
        give_active_elements_of_face(uint const & aFaceId) const;

        /**
         * Provides the position of an face
         *
         * @param[in] aFaceId            Face number.
         *
         * @param[out] Position         Position I,J,K.
         *
         */
        Mat<uint>
        give_face_position(uint & aFace)
        {
            return mBaseFace.give_face_position(aFace,mMeshData.ModelDim,mMeshData.NumberOfElementsPerDirection);
        }

        /**
         * Provides the ownership of an edge
         *
         * @param[in] aFaceId            Face number.
         *
         * @param[out] Rank        Proc rank
         *
         */
        uint
        give_face_owner(uint & aFaceId) const
        {
            return mBaseFace.give_face_owner(aFaceId,mMeshData.ModelDim,mMeshData.NumberOfElementsPerDirection,mProcData.DecompEdgeFace,mProcData.ProcNeighbor);
        }

        /**
         * Provides the a list, which share an face
         *
         * @param[in] aFaceId            Face number.
         *
         * @param[out] Share        Vector with procs, which share the face
         *
         */
        Mat<uint>
        give_face_share(uint & aFaceId) const
        {
            return mBaseFace.give_face_share(aFaceId,mMeshData.ModelDim,mMeshData.NumberOfElementsPerDirection,mProcData.DecompEdgeFace,mProcData.ProcNeighbor);
        }

        //Functionalities for B-spline basis and face----------------------------------
        /**
         * Provides the basis, which influence a face
         *
         * @param[in] aFaceId            Face number.

         *
         * @param[out] Nodes         Node id's, which are connected to a face
         *
         */
        Mat<uint>
        give_face_basis(uint const & aFaceId) const
        {
            return mFace.give_face_basis(aFaceId,mMeshData.ModelDim,mMeshData.NumberOfElementsPerDirection);
        }

        //Functionalities for Lagrange basis and face----------------------------------
        /**
         * Provides the basis, which influence a face
         *
         * @param[in] aFaceId            Face number.

         *
         * @param[out] Nodes         Node id's, which are connected to a face
         *
         */
        Mat<uint>
        give_face_lagrange_basis(uint const & aFaceId) const
        {
            return mLagrangeFace.give_face_basis(aFaceId,mMeshData.ModelDim,mMeshData.Polynomial,mMeshData.NumberOfElementsPerDirection);
        }

        /**
         * Provides the lagrange basis, which influence a face  (They are not really active, but always the lowest possible parent is provided )
         *
         * @param[in] aFaceId            Face number.

         *
         * @param[out] Nodes         Node id's, which are connected to a face
         *
         */
        Mat<uint>
        give_face_active_lagrange_basis(uint const & aFaceId) const;

        //Functions for TMatrix and IdField information-------------------------
        /**
         * Creates the elemental T-matrix for a 1,2 or 3 dimensional problem for p = 1,2,.... in a reordered fashion (classical FEM, counter clockwise numbering!)
         * Internal switch for non-truncated or truncated
         *
         * @param[in] aElementId                       Element Id.
         *
         * @param[out] aTMatrix        T-Matrix of the element will be saved internally in the struct
         * @param[out] aIdField        List of basis functions, which have support in this element, will be saved internally in the struct
         *
         */
        void
        give_Tmatrix_and_IdField_Reorder(
                uint const & aElementId);

        /**
         * Creates the elemental T-matrix for a 1,2 or 3 dimensional problem for p = 1,2,....
         * Internal switch for non-truncated or truncated
         *
         * @param[in] aElementId                       Element Id.
         *
         * @param[out] aTMatrix        T-Matrix of the element will be saved internally in the struct
         * @param[out] aIdField        List of basis functions, which have support in this element, will be saved internally in the struct
         *
         */
        void
        give_Tmatrix_and_IdField(
                uint const & aElementId);

        /**
         * Creates the elemental T-matrix for a 1,2 or 3 dimensional problem for p = 1,2,.... for the design field
         * Internal switch for non-truncated or truncated
         *
         * @param[in] aElementId                       Element Id.
         *
         * @param[out] aTMatrix        T-Matrix of the element will be saved internally in the struct
         * @param[out] aIdField        List of basis functions, which have support in this element, will be saved internally in the struct
         *
         */
        void
        give_Tmatrix_and_IdField_DesignVariables(
                uint const & aElementId);

        /**
         * Creates the elemental T-matrix for a 1,2 or 3 dimensional problem for p = 1,2,.... for the design field (This function nees to replace give_Tmatrix_and_IdField_DesignVariables
         * Internal switch for non-truncated or truncated
         *
         * @param[in] aElementId                       Element Id.
         *
         * @param[out] aTMatrix        T-Matrix of the element will be saved internally in the struct
         * @param[out] aIdField        List of basis functions, which have support in this element, will be saved internally in the struct
         *
         */
        void
        give_Tmatrix_and_IdField_DesignVariables_new(
                uint const & aElementId);

        /**
         * Creates the elemental T-matrix for a 1,2 or 3 dimensional problem for p = 1,2,.... for the design field of the last step
         * Internal switch for non-truncated or truncated
         *
         * @param[in] aElementId                       Element Id.
         *
         * @param[out] aTMatrix        T-Matrix of the element will be saved internally in the struct
         * @param[out] aIdField        List of basis functions, which have support in this element, will be saved internally in the struct
         *
         */
        void
        give_Tmatrix_and_IdField_DesignVariables_LastStep(
                uint const & aElementId);

        /**
         * Creates the nodal T-matrix for a 1,2 or 3 dimensional problem for p = 1,2,.... for a specific natural coordinate
         * Internal switch for non-truncated or truncated
         *
         * @param[in] aElementId                       Element Id.
         * @param[in] aNaturalCoordinate               A natural coordinate or a list of natural coordinates (rows = coordinates, cols = x, y, z)
         *
         * @param[out] aTMatrix        T-Matrix of the element will be saved internally in the struct
         * @param[out] aIdField        List of basis functions, which have support in this element, will be saved internally in the struct
         *
         */
        void
        give_Tmatrix_and_IdField_Specific_NaturalCoord(
                uint const & aElementId,
                Mat<real> const & aNaturalCoordinate);

        /**
         * Creates the projection matrix to map the solution from the Bsplines to the Lagrange nodes
         *
         * @param[in] aModelDim               Dimensions of the model
         * @param[in] aPolynomialLagrange      Polynomial degree of the lagrange basis
         * @param[in] aPolynomial              Polynomial degree of the Bspline basis
         *
         * @param[out] aTMatrix         Projection matrix
         *
         */
        Mat<real>
        give_projection_matrix_basis() const
        {
            return mTMatrix.give_projection_matrix_new( mMeshData.ModelDim,
                    mMeshData.PolynomialLagrange,
                    mMeshData.Polynomial );
        }

        /**
         * Creates the projection matrix to map the solution from the design Bsplines to the Lagrange nodes
         *
         * @param[in] aModelDim               Dimensions of the model
         * @param[in] aPolynomialLagrange      Polynomial degree of the lagrange basis
         * @param[in] aPolynomial              Polynomial degree of the Bspline basis
         *
         * @param[out] aTMatrix         Projection matrix
         *
         */
        Mat<real>
        give_projection_matrix_design_basis() const
        {
            return mTMatrix.give_projection_matrix_new( mMeshData.ModelDim,
                    mMeshData.PolynomialLagrange,
                    mMeshData.PolynomialDesign );
        }

        /**
         * Creates the nodal T-matrix for a 1,2 or 3 dimensional problem for p = 1,2,.... for a specific physical coordinate
         * Internal switch for non-truncated or truncated
         *
         * @param[in] aPhysicalCoordinate               A physical coordinate, which detects the active element and computes the natural coordinate
         *
         * @param[out] aTMatrix        T-Matrix of the element will be saved internally in the struct
         * @param[out] aIdField        List of basis functions, which have support in this element, will be saved internally in the struct
         *
         */
        void
        give_Tmatrix_and_IdField_Specific_PhysicalCoord(
                Mat<real> const & aPhysicalCoordinate);

        /**
         * Creates a vector with basis functions, which have support within an element
         * Internal switch for non-truncated or truncated
         *
         * @param[in] aElementId                       Element Id.
         *
         * @param[out] aIdField        List of basis functions, which have support in this element, will be saved internally in the struct
         *
         */
        Mat<uint>
        give_IdField(
                uint const & aElementId);

        /**
         * Creates a vector with design basis functions, which have support within an element
         * Internal switch for non-truncated or truncated
         *
         * @param[in] aElementId                       Element Id.
         *
         * @param[out] aIdField        List of design basis functions, which have support in this element, will be saved internally in the struct
         *
         */
        Mat<uint>
        give_Id_DesignField(
                uint const & aElementId);

        //Functions for Element and Basis Refinement----------------------------
        /**
         * Activates or refines elements, inherits the side sets, block sets and calculates the coordinates of the children
         *
         */
        void
        hierarchical_element_refinement();

        /**
         * Updates the elements in the outer layer to the passive state
         *
         *
         */
        void
        update_passive_outer_layer()
        {
            mRefinement.update_passive_outer_layer(
                    mMeshData.ModelDim,
                    mMeshData.Polynomial,
                    mMeshData.Level,
                    mMeshData.NumberOfElementsPerDirection,
                    mProcData.ElementRefinedListOnProc,
                    mElementData.ElementRefined );
        }

        /**
         * Activates or deactivates B-Spline basis functions
         *
         */
        void
        activate_basisfunction();

        /**
         *
         * Activates or deactivates Lagrange basis functions
         *
         */
        void
        activate_design_lagrange_basis();

        /**
         * Activates or deactivates basis functions (New functionality, which uses active,deactive,passive elements to active basis. Needs still to be proofed !!!!!)
         *
         */
        void
        activate_basisfunction_new();

        /**
         * Activates or deactivates design basis functions
         *
         *
         */
        void
        activate_basisdesignfunction();

        /** Find Elements, which have support with a basis function, which satisfy a nodal criteria, depending on a nodal field
         *
         * @param[in] aNodalField     Vector with nodal field data
         *
         */
        void
        refinement_criteria_nodal_field(
                map<uint, real> & aNodalField,
                const uint & aMaxRefinementLevel,
                const real & aLowerBound,
                const real & aUpperBound );

        /** Find Elements, which satisfy an elemental criteria, depending on an element field. The nodal field is needed to check if the element is solid or void
         *
         * @param[in] aElementField     Vector with elemental field data
         * @param[in] aNodalField     Vector with nodal field data
         *
         */
        void
        refinement_criteria_element_field(
                map<uint, real> & aElementField,
                map<uint, real> & aNodalField,
                const uint & aMaxRefinementLevel,
                const real & aLowerBound,
                const real & aUpperBound);
        /** Find Elements, which satisfy the signed distance field
         *
         * @param[in] aMaxLevelOfRefinement   Maximum level of refinement for the specific field
         * @param[in] aCandidateElementsToDeactivate    Vector with elements, which are a possible element which needs to be refined
         *
         */
        void
        find_elements_of_object(
                uint & aMaxLevelOfRefinement,
                const Mat<uint> & aCandidateElementsToDeactivate);

        /**
         * Find Elements in a stencil, which need to be refined
         *
         */
        void
        find_elements_in_stencil();

        /**
         * Creates a minimum refinement for the whole mesh
         *
         */
        void minimum_refinement()
        {
            if( mSettings.MinimumRefinement > 0)
            {
                mElementData.DeactivateElement = mRefinement.initial_refinement(
                        mMeshData.ModelDim,
                        mMeshData.NumberOfElementsPerDirection,
                        mSettings.MinimumRefinement,
                        mProcData.ElementListOnProcInit);
            }
        }

        /** Update basis values on the design elements
         *
         * @param[in] aNodalField     Nodal field on an equal or coarser mesh
         *
         * @param[out] aNodalField     Updated nodal field to the mesh of design elements
         *
         */
        void
        update_basis_values_on_design_elements(
                map<uint, real> & aNodalField);

        /** Update elemental values on the design elements. The elemental values are computed, based on nodal values
         *
         * @param[in] aNodalField     Nodal field on an equal or coarser mesh
         * @param[in] aElementField   Elemental field on an equal or coarser mesh
         *
         * @param[out] aElementField     Updated elemental field to the mesh of design elements
         *
         */
        void
        compute_element_field_on_design_elements(
                map<uint, real> & aNodalField,
                map<uint, real> & aElementField)
        {
            mRefinement.compute_element_field_on_design_elements(
                    mMeshData.ModelDim,
                    mMeshData.PolynomialDesign,
                    mMeshData.NumberOfElementsPerDirection,
                    mProcData.ElementListOnProcLastStepDesign,
                    mSettings.TruncatedBsplines,
                    aNodalField,
                    aElementField);
        }

        /**
         * Provides a list of active elements
         *
         * @param[out] ElementListOnProc      Updated list of elements, which are active on the current proc
         *
         */
        Mat<uint>
        give_element_on_proc()
        {
            mProcData.ElementListOnProc = mRefinement.give_active_elements(mMeshData.ModelDim,mMeshData.Polynomial,mMeshData.NumberOfElementsPerDirection,mMeshData.Level,mProcData.ElementListOnProcInit,mElementData.ElementActive);
            return mProcData.ElementListOnProc;
        }

        /**
         * Provides a list of active elements including an aura with a one element layer
         *
         * @param[out] ElementListOnProcWithAura      Updated list of elements, which are active on the current proc
         *
         */
        Mat<uint>
        give_element_on_proc_with_aura()
        {
            mProcData.ElementListOnProcWithAura = mRefinement.give_active_elements(mMeshData.ModelDim,mMeshData.Polynomial,mMeshData.NumberOfElementsPerDirection,mMeshData.Level,mProcData.ElementListOnProcWithAura,mElementData.ElementActive);
            return mProcData.ElementListOnProcWithAura;
        }

        /**
         * Creates a distributed nodeset ( a nodeset can be defined and the respective processor grabs the nodes which he owns)
         *
         * @param[in] aInNodeSet                    General nodeset which is user defined
         *
         * @param[out] aOutNodeSet   Distributed nodeset for each proc, which will be saved in the internal struct
         *
         */
        void
        set_nodeset(
                Mat<uint> & aNodeSet);

        /**
         * Creates a distributed side set. Each processor takes from the general side set the part he owns
         *

         * @param[in] aSideSet                      General sideset for the whole domain (Each processor grabs his necessary part)
         *
         * @param[out] aOutSideSet                  Distributed side set for each proc, which will be saved in the internal struct
         *
         */
        void
        set_sideset(
                Mat<uint> & aSideSet);

        /**
         * Updates the SideSet for the active elements/ basis functions
         *
         *
         */
        void
        update_sideset()
        {
            mMeshData.SideSetMTK = mSets.update_sideset(
                    mMeshData.ModelDim,
                    mMeshData.NumberOfElementsPerDirection,
                    mMeshData.Level,
                    mMeshData.NumSets,
                    mMeshData.SideSet,
                    mElementData.ElementActive);
        }
        /**
         * Creates a L2 projection for the design variables from the old mesh to the refined mesh of the design variables
         *
         * @param[in] aNodalField       A vector with a nodal field, which needs to be projected from the old to the new mesh
         *
         * @param[out] AbsDersVariables.dat    It creates an output file for femdoc
         *
         */
        void
        L2_projection(
                const map<uint, real> & aNodalField,
                const real aRhsScaling = 1.0,
                const real aRhsOffset  = 0.0 );

        /*void
        L2_projection_new(
            const map<uint, real> & aNodalField); */

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
                const uint & aDim,
                const uint & aPolynomial,
                const Mat<real> & aNodalField);

        static void
        RHS_for_L2_projection(
                const uint        & aDim,
                const uint        & aPolynomial,
                const Mat<real>   & aNodalField,
                      Mat< real > & aElementRHS );

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

        static void
        MassMatrix_for_L2_projection(
                const uint & aDim,
                const uint & aPolynomial,
                Mat<real>  & aElementMass );

        /**
         * Provides the active element, which has support for the specific coordinate
         *
         * @param[in] aCoordinate       A user specific coordinate
         *
         * @param[out] aElementId         Element Id of the active element
         *
         */
        uint
        give_active_element_for_coordinate(
                Mat<real> const & aCoordinate)
        {
            return mHMRElement.give_active_element_for_coordinate(mMeshData.ModelDim,mMeshData.NumberOfElementsPerDirection,mMeshData.Dimensions,mMeshData.DimensionsOffset,mSettings.PointOfOrigin,aCoordinate,mElementData.ElementActive);
        }

        /**
         * Provides the element on level zero, which has support for the specific coordinate (independent if it is active or deactive, needed for filter)
         *
         * @param[in] aCoordinate       A user specific coordinate
         *
         * @param[out] aElementId         Element Id of the active element
         *
         */
        uint
        give_element_for_coordinate_on_level_zero(
                Mat<real> & aCoordinate)
        {
            return mHMRElement.give_element_for_coordinate_on_level_zero(mMeshData.ModelDim,mMeshData.NumberOfElementsPerDirection,mMeshData.Dimensions,mMeshData.DimensionsOffset,mSettings.PointOfOrigin,aCoordinate);
        }

        /**
         * Provides basis functions in a vicinity of a given coordinate
         *
         * @param[in] aCoordinate       A user specific coordinate
         *
         * @param[out] aBasis           A list of basis functions, which have support in the element, which is in an active element
         *
         */
        Mat<uint>
        give_basis_vicinity_for_coordiante(Mat<real> & aCoordinate)
        {
            uint tElement = mHMRElement.give_active_element_for_coordinate(mMeshData.ModelDim,mMeshData.NumberOfElementsPerDirection,mMeshData.Dimensions,mMeshData.DimensionsOffset,mSettings.PointOfOrigin,aCoordinate,mElementData.ElementActive);
            return mHMRElement.give_basis_of_element(tElement,mMeshData.ModelDim,mMeshData.Polynomial,mMeshData.NumberOfElementsPerDirection);
        }

        /**
         * Provides a list of basis functions, which have support in a bounding box. The bounding box is represented with two coordinates
         *
         * @param[in] aCoordinate0       A user specific coordinate (smaller coordinate then aCoordinate1)
         * @param[in] aCoordinate1       A user specific coordinate (bigger coordinate then aCoordinate0)
         *
         * @param[out] aBasis           A list of basis functions, which have support in the bounding box
         *
         */
        Mat<uint>
        give_basis_vicinity_for_triangle(Mat<real> const & aCoordinate0,Mat<real> const & aCoordinate1);

        /**
         * Create all the data, which is necessary for the MTK
         *
         */
        //void
        //create_mesh_data();
        /**
         * Create all the data, which is necessary for the MTK
         *
         */
        void
        create_mesh_data_new();
        /**
         * Creates the MTK database to output a STK file
         *
         * @param[in] OutputFileName                Name of the output file
         *
         */
        void create_MTK_file(
                std::string const & aOutputFileName);

        /**
         * Creates the MTK database to output a STK file with field data
         *
         * @param[in] OutputFileName                Name of the output file
         * @param[in] aFieldData                    Field data for MTK for outputing elemental or nodal data
         * @param[in] aFieldName                    Name of the fiels
         * @param[in] aFieldRank                    Entity rank of the fields
         *
         */
        void
        create_MTK_file(
                std::string const & aOutputFileName,
                Cell<Mat<real>> & aFieldData,
                Cell<std::string> & aFieldName,
                Cell<enum EntityRank> & aFieldRank);

        /**
         * Creates the MTK database to output a STK file. It creates a mesh with data from the last continuation step
         *
         * @param[in] OutputFileName                Name of the output file
         *
         */
        void create_MTK_file_last_step(
                std::string const & aOutputFileName);

        /**
         * Creates the MTK database to output a STK file. It saves in addition to the standard output also the ADV as a nodal field
         *
         * @param[in] OutputFileName                Name of the output file
         *
         */
        void create_MTK_file_with_adv(
                std::string const & aOutputFileName);

        /**
         * Creates the MTK database to output a STK file. It saves in addition to the standard output also the ADV as a nodal field and the signed distance fields
         *
         * @param[in] OutputFileName                Name of the output file
         *
         */
        void create_MTK_file_with_adv_and_obj(
                std::string const & aOutputFileName );

        /**
         * Creates an output file for FemDoc with the elemental T-Matrix and Id-Field as ascii file
         *
         */
        void
        save_element_dof_connectivity_ascii();

        /**
         * Creates an output file for FemDoc with the elemental T-Matrix and Id-Field as binary file
         * as ascii file
         */
        void
        save_element_dof_connectivity_binary();

        /**
         * Creates an output file for FemDoc with the DOF variables and their interpolation per node
         *
         */
        void
        save_node_dof_connectivity_ascii();

        /**
         * Creates an output file for FemDoc with the DOF variables and their interpolation per node
         * as binary file
         */
        void
        save_node_dof_connectivity_binary();

        /**
         * Creates an output file for MORIS with the design variables and their interpolation per node
         *
         */
        void
        save_design_variables_for_moris_ascii();

        /**
         * Creates an output file for FEMDOC with the design variables and their interpolation per node
         *
         */
        void
        save_design_variables_for_femdoc_ascii();

        /**
         * Creates an output file for MORIS with the design variables and their interpolation per node
         *
         */
        void
        save_design_variables_for_moris_binary();
        /**
         * Creates an output file for FEMDOC with the design variables and their interpolation per node
         *
         */

        void
        save_design_variables_for_femdoc_binary();

        /**
         * Creates an output file for the continuation script to know in the next step the active design elements
         *
         */
        void
        save_active_design_elements();

        /**
         * Creates an output file for the continuation script to know in the next step the active FEM elements
         *
         */
        void
        save_active_fem_elements();

        /**
         * Creates an output file for the continuation script to know in the next step the active design basis functions
         *
         */
        void
        save_active_design_basis();

        /**
         * Creates an output file for the continuation script to know in the next step the active basis functions
         *
         */
        void
        save_active_basis();

        /**
         * Creates an output file for the continuation script to know in the next step the active basis functions
         *
         */
        void
        save_coordinate_list();

        /**
         * A wrapper to output all necessary files
         *
         */
        void
        save_all_data_files_for_continuation();
        /**
         * Reads a list of active FEM elements from an older step (e.g. continuation script)
         *
         */
        void
        read_active_fem_elements();

        /**
         * Reads a list of active FEM elements from an older step (e.g. continuation script)
         *
         */
        void
        read_active_design_elements();

        /**
         * Reads a list of active design basis functions from an older step (e.g. continuation script)
         *
         */
        void
        read_active_design_basis();

        /**
         * Reads a list of active basis functions from an older step (e.g. continuation script)
         *
         */
        void
        read_active_basis();

        /**
         * Reads the nodal ID design basis functions and T-Matrix from the designvariables.data from an older step ( e.g. continuation script)
         *
         */
        void
        read_design_variables_ascii();

        /**
         * Reads a list of coordinate numbers of the last step
         *
         */
        void
        read_coordinate_list();

        /**
         * Reads data from the AbsDesVariables file and puts user specific calculations on the different fields
         *
         */
        void
        read_absdesvariables_file();

        /**
         * A wrapper to read all necessary files
         *
         */
        void
        read_all_data_files_for_continuation();
        /**
         * All bitsets are send to proc 0, which compares data and save on tBitsetAlpha only those, where both vectors have at the same position a one. Afterwards this tBitsetAlpha is broadcasted to all procs
         *
         * @param[in] aBitset       Bitset with active elements or basis functions
         *
         * @param[out] aBitset      Bitset with active elements or basis functions after comparing with all procs
         *
         */
        void
        read_sdf_for_active_design_basis();

        void
        broadcast_bitset_logical_and(
                BoostBitset & aBitset);

        /**
         * All bitsets are send to proc 0, which compares data and save on tBitsetAlpha all bits with a one. Afterwards this tBitsetAlpha is broadcasted to all procs
         *
         * @param[in] aBitset       Bitset with active elements or basis functions
         *
         * @param[out] aBitset      Bitset with active elements or basis functions after comparing with all procs
         *
         */
        void
        broadcast_bitset_logical_or(
                BoostBitset & aBitset);

        /**
         * Gather a message to proc 0, take the max of the proc values and broadcast the max.
         *
         * @param[in] aElementList       aElementList.Level
         *
         * @param[out] aElementList      aElementList.Level (updated)
         *
         */
        void
        gather_value_and_bcast_max(
                uint & aMessage)
        {
            mMPI.gather_value_and_bcast_max( aMessage );
        }

        //Floodfill algorithms
        /**
         * Create a nodal vector with plus or minus
         *
         * @param[in] aNodalField        A vector with plus or minus
         *
         * @param[out] aBasisList       Updates the basis functions: active/deactive, node sets
         * @param[in] aElementList      Updates the elements: (see top: ElementList_struc)
         *
         */
        void floodfill_for_STL(Mat<real> & aNodalField);

        // Functions to have access to struc data-------------------------------
        /**
         * Function to return the element dimensions
         */
        Mat<real>
        give_element_dimensions()
        {
            return mMeshData.ElementLength;
        }

        /**
         * Function to return the number of elements in each direction
         *
         * @TODO consider renaming this to give_number_elements_per_direction
         *
         */
        Mat<uint>
        give_set_number_elements_per_direction()
        {
            return mMeshData.NumberOfElementsPerDirection;
        }

        /**
         *  Function to return the model dimension
         */
        uint
        give_modeldim() const
        {
            return mMeshData.ModelDim;
        }

        /**
         *  Function to return the polynomial
         */
        uint
        give_polynomial() const
        {
            return mMeshData.Polynomial;
        }

        /**
         *  Function to return the domain dimensions (domain range)
         */
        Mat<real>
        give_domain_dimensions()
        {
            return mMeshData.Dimensions;
        }

        /**
         * Function to return the point of origin of the domain
         */
        Mat<real>
        give_point_of_origin()
        {
            return mSettings.PointOfOrigin;
        }

        /**
         * Function to return the prescribed side sets
         */
        Mat<uint>
        give_side_set()
        {
            return mMeshData.SideSet;
        }

        /**
         * Function to return the domain, which needs to be deactive from the beginning
         * (returns a 2x3 maxtrix, first point is the lower point of the bounding box and the second the upper point)
         */
        Mat<real>
        give_deactive_domain()
        {
            return mMeshData.DeleteDomain;
        }

        /**
         * Function which returns the number of the buffer layer
         */
        /*uint
        give_buffer(){
            return mMeshData.Buffer;
        }*/

        /**
         * Function to set the point of origin
         *
         * @param[in] aOrigin Mat<real> containing the coordinates of first
         *                              node of the calculation domain
         */
        void
        set_point_of_origin(
                const Mat<real>& aOrigin)
        {
            mSettings.PointOfOrigin = aOrigin;
        }

        /**
         * Function to set the numbers of elements per direction
         *
         * @param[in] aNumberOfElementsPerDirection
         */
        void
        set_number_elements_per_direction(
                Mat<uint> & aNumberOfElementsPerDirection)
        {
            mMeshData.NumberOfElementsPerDirection = aNumberOfElementsPerDirection;
        }

        /**
         * Function to set the domain dimensions (domain range)
         *
         * @param[in] aDimensions width, height and depth of calculation domain
         */
        void
        set_domain_dimensions(
                Mat<real> & aDimensions)
        {
            mMeshData.Dimensions = aDimensions;
        }

        /**
         * Function to set the deactive domain from the beginning
         * First point is the lower and second point the upper point of a bounding box
         */
        void
        set_deactive_domain(
                Mat<real> & aRemovePoint1,
                Mat<real> & aRemovePoint2 )
        {
            mMeshData.DeleteDomain.row(0) = aRemovePoint1.row(0);
            mMeshData.DeleteDomain.row(1) = aRemovePoint2.row(0);
        }

        /**
         * Function to set the deactive domain from the beginning (input is a 2x3 matrix)
         */
        void
        set_deactive_domain(
                Mat<real> & aRemovePoint )
        {
            mMeshData.DeleteDomain = aRemovePoint;
        }

        /**
         * Function to set the side set
         */
        void
        set_side_set(
                Mat<uint> & aSideSet)
        {
            mMeshData.SideSet = aSideSet;
        }

        /**
         * Function to set the point of origin
         *
         * @param[in] aOrigin Mat<real> containing the coordinates of first
         *                              node of the calculation domain
         */
        void
        set_timing_tests()
        {
            mSettings.TimingTest = true;
        }

        /**
         * Copy data from a dummy mesh to the used mesh from a parser
         */
        void
        copy_parser_information(
                Hierarchical_Mesh_Main& aHierarchicalMeshIn,
                Hierarchical_Mesh_Main& aHierarchicalMeshOut);

        /**
         * Provides a list of langrange entities, which are owned and shared by the current proc
         *
         * @param[in] aEntityRank        Rank of the entity
         *
         * @param[out] tEntities         Entities, which are owned and share by the current processor
         */
        Mat<uint>
        get_entities_owned_and_shared_by_current_proc(
                enum EntityRank   aEntityRank) const;

        /**
         * Provides a list of bspline, which are owned and shared by the current proc
         *
         * @param[in] aEntityRank        Rank of the entity
         *
         * @param[out] tEntities         Entities, which are owned and share by the current processor
         */
        Mat<uint>
        get_bspline_entities_owned_and_shared_by_current_proc(
                enum EntityRank   aEntityRank) const;

        /**
         * Provides the number of entities, which are owned and shared by the current proc
         *
         * @param[in] aEntityRank        Rank of the entity
         *
         * @param[out] tNumEntities      Number of entities, which are owned and shared by the current processor
         */
        uint
        get_num_entities_owned_and_shared_by_current_proc(
                enum EntityRank   aEntityRank) const;

        /**
         * Provides a list of entities, which are owned and in the aura of the current proc
         *
         * @param[in] aEntityRank        Rank of the entity
         *
         * @param[out] tEntities      Entities, which are owned and in the aura of the current processor
         */
        Mat<uint>
        get_entities_universal(
                enum EntityRank   aEntityRank) const;

        /**
         * Provides the number of entities, which are owned and in the aura of the current proc
         *
         * @param[in] aEntityRank        Rank of the entity
         *
         * @param[out] tNumEntities      Number of entities, which are owned and in the aura of the current processor
         */
        uint
        get_num_entities_universal(
                enum EntityRank   aEntityRank) const;

        /**
         * Provides a list of entities, which are owned by the current proc
         *
         * @param[in] aEntityRank        Rank of the entity
         *
//         * @param[out] tEntities      Entities, which are owned by the current processor
         */
        Mat<uint>
        get_entities_owned_current_proc(
                enum EntityRank   aEntityRank) const;

        /**
         * Provides the number of entities, which are owned by the current proc
         *
         * @param[in] aEntityRank        Rank of the entity
         *
         * @param[out] tNumEntities      Number of entities, which are in the owned by the current processor
         */
        uint
        get_num_entities_owned_current_proc(
                enum EntityRank   aEntityRank) const;

        /**
         * Provides a list of entities, which are in the aura of the current proc
         *
         * @param[in] aEntityRank        Rank of the entity
         *
         * @param[out] tEntityAura       Entities, which are in the aura of the current processor
         */
        Mat<uint>
        get_entities_in_aura(
                enum EntityRank   aEntityRank) const;

        /**
         * Provides the number of entities, which are in the aura of the current processor
         *
         * @param[in] aEntityRank        Rank of the entity
         *
         * @param[out] tNumEntities      Number of entities, which are in the aura of the current processor
         */
        uint
        get_num_entities_in_aura(
                enum EntityRank   aEntityRank) const;

        /**
         * Provides the owner of a given entity
         *
         * @param[in] aEntityID          Id of the entity
         * @param[in] aEntityRank        Rank of the entity
         *
         * @param[out] tProcRank        Rank of the processor, who owns the entity
         *
         */
        uint
        parallel_owner_rank_by_entity_id(
                uint aEntityID,
                enum EntityRank aEntityRank ) const;

        /**
         * Provides a list of processors, who share a specific entity
         *
         * @param[in] aEntityID          Id of the entity
         * @param[in] aEntityRank        Rank of the entity
         *
         * @param[out] tShareProcs       Rank of processors, who share the entity
         *
         */
        Mat<uint>
        get_procs_sharing_entity_by_id(
                uint aEntityID,
                enum EntityRank aEntityRank ) const;

        /**
         * Provides a list of processors, who share a specific entity
         *
         * @param[in] aGFile          Name of exo file
         *
         * @param[out] aLevelsetfield     Levelset field from data
         * @param[out] aDensityfield      Density field from data
         *
         */
        void
        read_data_from_exo_file(
                std::string & aGFile);

        /**
         * Provides a vector for reordering (needed for paraview)
         *
         */
        Mat<uint>
        give_vector_for_reorder() const;
// -------------------------------------------------------------------------------------------------
// functions for new basis numbering
// -------------------------------------------------------------------------------------------------

        uint convert_lagrange_id_to_bspline_id(
                const uint aMaxLevel,
                const uint aGlobalLagrangeID) const ;

// -------------------------------------------------------------------------------------------------

        uint convert_lagrange_id_to_design_bspline_id(
                const uint aMaxLevel,
                const uint aGlobalLagrangeID) const ;

// -------------------------------------------------------------------------------------------------
        /**
        * generate maps needed for mesh data calcuation and output
        */
        void
        create_maps();

// -------------------------------------------------------------------------------------------------
// element and basis consistency checks
// -------------------------------------------------------------------------------------------------
        /**
         * make sure that each position in them mesh has not more than one active B-Spline basis
         */
        void
        check_active_basis();

// -------------------------------------------------------------------------------------------------
        /**
         * make sure that if an element is active, its parent is not active
         */
        void
        check_active_elements();

// -------------------------------------------------------------------------------------------------
// various functions needed for SDF
// -------------------------------------------------------------------------------------------------
        /**
         * calls L2 projection to map ADV on new field
         */
        void
        calculate_ADV_on_new_field();

// -------------------------------------------------------------------------------------------------

        /**
         * converts the data generated by SDF so that refinement_criteria_element_field and
         * refinement_criteria_nodal_field can understand
         */
        void
        create_field_maps_from_sdf(
                const Mat< real > & aSignOfSDF ,
                const Mat< uint > & aElementsAtSurface,
                const Mat< uint > & aElementsInVolume );

// -------------------------------------------------------------------------------------------------

        void
        create_nodal_fields_from_sdf();

// -------------------------------------------------------------------------------------------------

        /*
         * checkes if a SDF is to be used for refinement
         */
        void
        init_sdf_bitset();

// -------------------------------------------------------------------------------------------------
// various new functions
// -------------------------------------------------------------------------------------------------

        void
        check_and_init_symmetry();

// -------------------------------------------------------------------------------------------------

        uint
        give_symmetry_master_of_basis( const uint& aBasisID ) const;

// -------------------------------------------------------------------------------------------------

        uint
        give_symmetry_master_of_element( const uint& aElementID ) const;

// -------------------------------------------------------------------------------------------------

        uint
        give_symmetry_slave_of_element( const uint& aElementID ) const;

// -------------------------------------------------------------------------------------------------

        bool
        basis_is_symmetry_master( const uint& aBasisID ) const;

// -------------------------------------------------------------------------------------------------

        /**
         * makes sure that symmetry is maintained when mesh is refined
         */
        void
        make_deactivate_list_symmetric();

// -------------------------------------------------------------------------------------------------

        /**
         * makes sure that the connectivities in output for FEMDOC are consistent
         * exits MORIS with an error if not.
         */
        void
        check_dof_consistency();

// -------------------------------------------------------------------------------------------------

        /**
         * needed to make SDF symmetric
         */
        void
        activate_symmetry_design_bitset();

// -------------------------------------------------------------------------------------------------

    };
}

#endif /* SRC_MESH_CL_HIERARCHICAL_MESH_MAIN_HPP_ */

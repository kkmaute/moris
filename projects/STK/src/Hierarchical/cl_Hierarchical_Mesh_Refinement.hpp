/*
 * cl_Hierarchical_Mesh_Refinement.hpp
 *
 *  Created on: Dec 21, 2017
 *      Author: gleim
 */

#ifndef SRC_MESH_CL_HIERARCHICAL_MESH_REFINEMENT_HPP_
#define SRC_MESH_CL_HIERARCHICAL_MESH_REFINEMENT_HPP_

#include "linalg.hpp"
#include "cl_Base_Mat.hpp" // LNA/src
#include "cl_Mat.hpp" // LNA/src
#include "cl_Map.hpp" // CON/src
#include "cl_BoostBitset.hpp" // CON/src
#include "fn_unique.hpp" // LNA/src
#include "fn_histc.hpp" // LNA/src
#include "fn_sum.hpp" // LNA/src
#include "cl_Bspline.hpp" // MOD/src
#include "cl_Base_Mesh_Element.hpp"
#include "cl_Hierarchical_Mesh_Element.hpp"
#include "cl_Hierarchical_Mesh_Basis.hpp"
#include "cl_Hierarchical_Mesh_MPI.hpp"
#include "cl_Hierarchical_Mesh_TMatrix.hpp"

namespace moris
{

    class Hierarchical_Mesh_Refinement
    {
    protected:

    public:
        //Create Object of MPI
        Hierarchical_Mesh_MPI mMPI;
        //Create Object of TMatrix
        Hierarchical_Mesh_TMatrix mTMatrix;
        //Create Object of Basis
        Hierarchical_Mesh_Basis mBasis;
        //Create Object of Element
        Hierarchical_Mesh_Element mHMRElement;
        //Create Object of base Element
        Base_Mesh_Element mBaseElement;
        /**
         * Hierarchical_Mesh constructor
         */
        Hierarchical_Mesh_Refinement()
    {
    }

        /**
         * Hierarchical_Mesh destructor.
         */
        ~Hierarchical_Mesh_Refinement() = default;

        /**
         * Activates or deactivates elements, inheret the side sets, block sets and calculates the coordiantes of the children
         *
         * @param[in] aModelDim                Number of dimensions.
         * @param[in] aPolynomial          Polynomial degree of the basis functions.
         * @param[in] aNumberOfElementsPerDirection        NumElements=[Number of elements in x,y,z direction] 1x3.
         * @param[in] aDeactivateElement        Vector of all elements, which need to be deactivated.
         * @param[in] aElementActive      elements: active and deactive elements in bitset
         *
         * @param[out] aElementActive       Updates the elements: active/deactive in bitset
         *
         */
        void hierarchical_element_refinement(
                uint const & aModelDim,
                uint const & aPolynomial,
                Mat<uint> const & aNumberOfElementsPerDirection,
                Mat<uint> & aDeactivateElement,
                BoostBitset & aElementActive,
                BoostBitset & aElementRefined);

        /**
         * Flags the elements in the outer layer to be passive
         *
         * @param[in] aModelDim                          Number of dimensions.
         * @param[in] aPolynomial                        Polynomial degree of the basis functions.
         * @param[in] aLevel                             finest level
         * @param[in] aNumberOfElementsPerDirection        NumElements=[Number of elements in x,y,z direction] 1x3.
         * @param[in] aElementRefinedListOnProc         Vector of all elements, which are in the outer layer on the coarsest mesh
         * @param[in] aElementRefined                   Bitset with passive elements
         *
         * @param[out] aElementRefined                 Updated bitset with passive elements in outer layer
         *
         */
        void
        update_passive_outer_layer(
                uint const & aModelDim,
                uint const & aPolynomial,
                uint const & aLevel,
                Mat<uint> const & aNumberOfElementsPerDirection,
                Mat<uint> const & aElementRefinedListOnProc,
                BoostBitset & aElementRefined);

        /**
         * Activates or deactivates basis functions
         *
         * @param[in] aModelDim                Number of dimensions.
         * @param[in] aPolynomial          Polynomial degree of the basis functions.
         * @param[in] aNumberOfElementsPerDirection        NumElements=[Number of elements in x,y,z direction] 1x3.
         * @param[in] aBasisActive        Basis: active and deactive elements in bitset
         * @param[in] aElementActive      Elements: active and deactive elements in bitset
         * @param[in] ElementListOnProc   List of elements, which are active on the current proc
         *
         * @param[out] aBasisActive       Updates the basis function: active/deactive in bitset
         *
         * Exception:
         *       * elements on the edge outside of the domain are always refined
         */
        void activate_basisfunction(
                uint const & aModelDim,
                uint const & aPolynomial,
                Mat<uint> const & aNumberOfElementsPerDirection,
                BoostBitset & aBasisActive,
                BoostBitset const & aElementActive,
                Mat<uint> const & aElementListOnProc,
                Mat<uint> & aActiveBasisList);

        /**
         * Activates or deactivates basis functions, inherit the node sets
         *
         * @param[in] aModelDim                Number of dimensions.
         * @param[in] aPolynomial          Polynomial degree of the basis functions.
         * @param[in] aNumberOfElementsPerDirection        NumElements=[Number of elements in x,y,z direction] 1x3.
         * @param[in] aBasisActive        Basis: active and deactive elements in bitset
         * @param[in] aElementActive      Elements: active and deactive elements in bitset
         * @param[in] ElementListOnProc   List of elements, which are active on the current proc
         *
         * @param[out] aBasisActive       Updates the basis function: active/deactive in bitset
         *
         * A basis function is
         *       * Active (aka awake), if all its elements are either active or refined, and at least one is active
         *       * Refined (aka resting aka passive), if all its elements are refined
         *       * Deactive (aka dead aka unborn), otherwise
         *
         */
        void
        activate_basisfunction_new(
                uint const & aModelDim,
                uint const & aPolynomial,
                Mat<uint> const & aNumberOfElementsPerDirection,
                uint const & aLevel,
                BoostBitset & aBasisActive,
                BoostBitset & aBasisRefined,
                BoostBitset const & aElementActive,
                BoostBitset const & aElementRefined,
                Mat<uint> const & aElementListOnProc,
                Mat<uint> const & aElementRefinedListOnProc,
                Mat<uint> & aActiveBasisList);

        /**
         * Provides a list of active elements
         *
         * @param[in] aModelDim                    Number of dimensions.
         * @param[in] aPolynomial             Polynomial degree of the basis functions.
         * @param[in] aNumberOfElementsPerDirection            NumElements=[Number of elements in x,y,z direction] 1x3.
         * @param[in] aLevel                  Current max Level
         * @param[in] ElementListOnProcInit   List of elements, which are active on the initial mesh
         * @param[in] aElementActive          Elements: active and deactive elements in bitset
         *
         * @param[out] ElementListOnProc      Updated list of elements, which are active on the current proc
         *
         */
        Mat<uint>
        give_active_elements(
                uint const & aModelDim,
                uint const & aPolynomial,
                Mat<uint> const & aNumberOfElementsPerDirection,
                uint const & aLevel,
                Mat<uint> const & aElementListOnProcInit,
                BoostBitset const & aElementActive) const;

        /**
         * Provides a list of deactivated elements
         *
         * @param[in] aModelDim                Number of dimensions.
         * @param[in] aPolynomial          Polynomial degree of the basis functions.
         * @param[in] aNumberOfElementsPerDirection        NumElements=[Number of elements in x,y,z direction] 1x3.
         * @param[in] aElementActive      Sparsity matrix with list of active/deactive elements.
         * @param[in] aBasisActive        Sparsity matrix with list of active/deactive basis functions.
         * @param[in] aDeactivateElement        Vector of all elements, which need to be deactivated.
         *
         * @param[out] aDeactivateElement       Updated list of deactivated elements (linear dependencies, ...) .
         *
         */
        void
        give_deactivated_elements(
                uint & aModelDim,
                uint & aPolynomial,
                Mat<uint> & aNumberOfElementsPerDirection,
                Mat<uint> & aDeactivateElements,
                BoostBitset & aElementActive,
                Mat<uint> & aElementListOnProcInit);

        /**
         * Find Elements in a stencil, which need to be refined
         *
         * @param[in] aModelDim                Number of dimensions.
         * @param[in] aPolynomial          Polynomial degree of the basis functions.
         * @param[in] aNumberOfElementsPerDirection        NumElements=[Number of elements in x,y,z direction] 1x3.
         * @param[in] aElementListOnProc   A list of elements, which are active
         * @param[in] aElementField        A vector which an elemental property, which is the criteria
         * @param[in] aFeatureLowerBound   A threshold, which decides wheter a element needs to be deactivated or not
         * @param[in] aFeatureResolution   An element gets deactivated, if this threshold is valid
         *
         * @param[out] aElementList       A list of elements, which need to be deactivated
         *
         */
        Mat<uint>
        find_elements_in_stencil(
                uint & aModelDim,
                uint & aPolynomial,
                Mat<uint> & aNumberOfElementsPerDirection,
                BoostBitset & aElementActive,
                Mat<uint> & aElementListOnProc,
                map<uint, real> & aElementField,
                real & aFeatureLowerBound,
                uint & aFeatureResolution,
                const uint      & aBufferElements,
                const bool      & aAdaptiveBufferLayer,
                const bool      & aStaircasebuffer);

        /* Find Elements, which have support with a basis function above a threshold, which need to be refined. And the second criteria are intersected Elements ( Level set zero contour)
         *
         * @param[in] aElementList     Data of ElementList
         * @param[in] aBasisList       Data of BasisList
         *
         * @param[out] aElementList       ElementList.DeactivateElements
         *
         */
        void
        find_elements_for_refinement(
                uint & aModelDim,
                uint & aPolynomialDesign,
                Mat<uint> & aNumberOfElementsPerDirection,
                Mat<uint> & aElementListLastStepDesign,
                Mat<uint> & aDeactivateElements,
                Mat<real> & aElementField,
                BoostBitset & aElementActive,
                BoostBitset & aBasisActive,
                Mat<real> & aNodalField,
                bool & aTruncatedBsplines,
                real & aDensityLowerBound,
                real & aDensityIncrementBound,
                real & aDensityUpperBound,
                uint & aMaxDesignLevelOfRefinement,
                bool & aLvlSetMethod,
                uint & aMaxLevelSetLevelOfRefinement,
                uint & aBufferElements,
                bool & aAdaptiveBufferLayer,
                bool & aStaircasebuffer);

        /* Find Elements, which have support with a basis function above a threshold, which need to be refined. And the second criteria are intersected Elements ( Level set zero contour)
         *
         * @param[in] aElementList     Data of ElementList
         * @param[in] aBasisList       Data of BasisList
         *
         * @param[out] aElementList       ElementList.DeactivateElements
         *
         */
        void
        refinement_criteria_nodal_field(
                uint & aModelDim,
                uint & aPolynomialDesign,
                Mat<uint> & aNumberOfElementsPerDirection,
                Mat<uint> & aElementListLastStepDesign,
                Mat<uint> & aDeactivateElements,
                BoostBitset & aElementActive,
                BoostBitset & aBasisActive,
                map<uint, real> & aNodalField,
                const bool & aTruncatedBsplines,
                const real & aNodalLowerBound,
                const real & aNodalUpperBound,
                const uint & aMaxDesignLevelOfRefinement,
                const uint & aBufferElements,
                const bool & aAdaptiveBufferLayer,
                const bool & aStaircasebuffer);

        /* Find Elements, which have support with a basis function above a threshold, which need to be refined. And the second criteria are intersected Elements ( Level set zero contour)
          *
          * @param[in] aElementList     Data of ElementList
          * @param[in] aBasisList       Data of BasisList
          *
          * @param[out] aElementList       ElementList.DeactivateElements
          *
          */
         void
         refinement_criteria_element_field(
                 uint & aModelDim,
                 uint & aPolynomialDesign,
                 Mat<uint> & aNumberOfElementsPerDirection,
                 Mat<uint> & aElementListLastStepDesign,
                 Mat<uint> & aDeactivateElements,
                 map<uint, real> & aElementField,
                 BoostBitset & aElementActive,
                 BoostBitset & aBasisActive,
                 map<uint, real> & aNodalField,
                 const bool & aTruncatedBsplines,
                 const real & aDensityLowerBound,
                 const real & aDensityIncrementBound,
                 const real & aDensityUpperBound,
                 const uint & aMaxDesignLevelOfRefinement,
                 const uint & aBufferElements,
                 const bool & aAdaptiveBufferLayer,
                 const bool & aStaircasebuffer);

        /* Find the intersected elements of an object file
         *
         * @param[in] aElementList     Data of ElementList
         * @param[in] aBasisList       Data of BasisList
         *
         * @param[out] aElementList       ElementList.DeactivateElements
         *
         */

        // call me from Mesh Main : line 737
        void
        find_elements_of_object(
                const uint & aModelDim,
                const uint & aPolynomial,
                const uint & aMaxLevelOfRefinement,
                const Mat<uint> & aNumberOfElementsPerDirection,
                BoostBitset & aElementActive,
                const uint & aBufferElements,
                const uint & aBufferLevel,
                bool_t & aAdaptiveBufferLayer,
                bool_t & aStaircasebuffer,
                const Mat<uint> & aCandidateElementsToDeactivate,
                Mat<uint>& aElementsToDeactivate);
        /* Create an inital refinement for the first step
         *
         * @param[in] aElementList     Data of ElementList
         * @param[in] aBasisList       Data of BasisList
         *
         * @param[out] aElementList       ElementList.DeactivateElements
         *
         */
        Mat<uint>
        initial_refinement(
                uint const & aModelDim,
                Mat<uint> const & aNumberOfElementsPerDirection,
                uint const & aInitialRefinement,
                Mat<uint> & aElementListOnProcInit);

        /* Creates a buffer layer for a list of elements, which need to be deactivated
         *
         * @param[in] aElementList     Data of ElementList
         * @param[in] aBasisList       Data of BasisList
         *
         * @param[out] aElementList       ElementList.DeactivateElements
         *
         */
        void
        create_buffer_layer(
                const uint & aModelDim,
                const Mat<uint> & aNumberOfElementsPerDirection,
                const uint & aLevel,
                const uint & aBufferElements,
                const bool & aAdaptiveBufferLayer,
                const bool & aStaircasebuffer,
                Mat<uint> & aDeactivatElement,
                BoostBitset & aElementActive);

        /* Searches for parent elements, which also need to be deactivated
         *
         * @param[in] aElementList     Data of ElementList
         * @param[in] aBasisList       Data of BasisList
         *
         * @param[out] aElementList       ElementList.DeactivateElements
         *
         */
        void
        add_parents_to_list(
                const uint & aModelDim,
                const Mat<uint> & aNumberOfElementsPerDirection,
                Mat<uint> & aRefineElement,
                BoostBitset & aElementActiveDummy);

        /* Update basis values on the design elements
         *
         * @param[in] aElementList     Data of ElementList
         * @param[in] aBasisList       Data of BasisList
         *
         * @param[out] aElementList       ElementList.DeactivateElements
         *
         */
        void
        update_basis_values_on_design_elements(
                uint & aModelDim,
                uint & aPolynomialDesign,
                Mat<uint> & aNumberOfElementsPerDirection,
                Mat<uint> & aElemLocaltoGlobalLastStepDesign,
                bool & aTruncatedBsplines,
                BoostBitset & aDesignBSplineActiveLastStep,
                map<uint, real> & aNodalField,
                Cell<Mat<real>> & aTMatrixParentChild);

        /* Update elemental value with the solution of basis values
         *
         * @param[in] aElementList     Data of ElementList
         * @param[in] aBasisList       Data of BasisList
         *
         * @param[out] aElementList       ElementList.DeactivateElements
         *
         */
        void
        compute_element_field_on_design_elements(
                uint & aModelDim,
                uint & aPolynomialDesign,
                Mat<uint> & aNumberOfElementsPerDirection,
                Mat<uint> & aElemLocaltoGlobalLastStepDesign,
                bool & aTruncatedBsplines,
                map<uint, real> & aNodalField,
                map<uint, real> & aElementField);

        // this function is not needed at the moment.
        /* void
        finish_refinement_list(
                Mat< uint >           & aDeactivateElements,
                const BoostBitset     & aElementActive,
                const uint            & aModelDim,
                const uint            & aPolynomialDesign,
                const Mat<uint>       & aNumberOfElementsPerDirection,
                const uint            & aBufferElements,
                const bool            & aAdaptiveBufferLayer,
                const bool            & aStaircasebuffer); */
    };

}

#endif /* SRC_MESH_CL_HIERARCHICAL_MESH_REFINEMENT_HPP_ */

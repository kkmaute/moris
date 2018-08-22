/*
 * cl_Lagrange_Filter.hpp
 *
 *  Created on: Mar 7, 2018
 *      Author: gleim
 */

#ifndef SRC_MESH_CL_LAGRANGE_FILTER_HPP_
#define SRC_MESH_CL_LAGRANGE_FILTER_HPP_

#include "algorithms.hpp"
#include "linalg.hpp"
#include "chronos.hpp"
#include "cl_Base_Mat.hpp" // LNA/src
#include "cl_Mat.hpp" // LNA/src
#include "cl_Cell.hpp" // CON/src
#include "cl_Map.hpp" // CON/src
#include "fn_unique.hpp" // LNA/src
#include "fn_find_unique.hpp" // LNA/src
#include "fn_find.hpp" // LNA/src
#include "fn_sum.hpp" // LNA/src
#include "fn_norm.hpp"
#include "op_elemwise_mult.hpp" // LNA/src
#include "cl_BoostBitset.hpp" // CON/src
#include "cl_Hierarchical_Mesh_Element.hpp"
#include "cl_Lagrange_Basis.hpp"
#include "cl_Lagrange_Element.hpp"
#include "cl_Base_Mesh_Element.hpp"

namespace moris
{

    uint gCandidateCounter = 0;
    uint gCandidateHitCounter = 0;

    class Lagrange_Filter
    {
    protected:

    public:
        //Create Object of Basis
        Lagrange_Basis mLagrangeBasis;
        //Create Object of Element
        Lagrange_Element mLagrangeElement;
        //Create Object of Element
        Hierarchical_Mesh_Element mHMRElement;
        //Create Object of Element
        Base_Mesh_Element mBaseElement;

        /**
         * Hierarchical_Mesh constructor
         */
        Lagrange_Filter()
    {
    }

        /**
         * Hierarchical_Mesh destructor.
         */
        ~Lagrange_Filter() = default;

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
                uint const & aDim,
                uint const & aPolynomialDesign,
                Mat<uint> const & aNumElements,
                uint const & aLevel,
                real const & aFilterRadius,
               BoostBitset & aElementActiveDesign,
               const BoostBitset & aDesignBasisActive,
               const BoostBitset & aLagrangeBasisActive,
               const uint        & aNumberOfActiveDesignLagrangeBasis,
               //const map<uint, uint> & aDesignBSplineListMap,
               Mat<real> const & aDimensions,
               Mat<real> const & aDimensions_Offset,
               Mat<real> const & aDimensionsOriginal,
               Mat<uint> & aIdFieldDesignNodalField,
               Mat<real> & aTMatrixDesignNodalField,
               bool const & aFilterLevelWeight,
               bool const & aPerformNormalization,
               Mat<real> const & aPointOfOrigin,
               Mat<uint> & IdFieldDesignNodalList);

//        /**
//         * Provides a list of basis functions, which are in a specified circle
//         *
//         * @param[in] aBasis             Basis function, which is used to search for neighbours
//         * @param[in] aElementList       Struc which contains the information about the mesh
//         * @param[in] aBasisList       Struc which contains the information about the basis functions
//         *
//         * @param[out] Mat             Output is a matrix with three columns. First column are a list of elements, second and third column are the distance and the ratio level
//         *
//         */
//        Mat<real>
//        Filter_for_smoothing(
//                uint const & aBasis,
//                uint const & aDim,
//                uint const & aPolynomialDesign,
//                Mat<uint> const & aNumElements,
//                uint const & aLevel,
//                real const & aFilterRadius,
//               BoostBitset & aElementActiveDesign,
//               BoostBitset & aDesignBSplineActive,
//               Mat<real> const & aDimensions,
//               Mat<real> const & aDimensions_Offset);

        /**
         * Provides a list of basis functions, which are in a specified circle or sphere
         *
         * @param[in] aBasis           Basis function, which is used to search for neighbors
         * @param[in] aElementList     Struc which contains the information about the mesh
         * @param[in] aBasisList       Struc which contains the information about the basis functions
         *
         * @param[out] Mat             Output is a matrix with three columns. First column are a list of elements, second and third column are the distance and the ratio level
         *
         */
        void
        filter_for_smoothing(
                const uint               & aBasis,
                const uint               & aDim,
                const uint               & aPolynomialDesign,
                const Mat<uint>          & aNumElements,
                const uint               & aLevel,
                const real               & aFilterRadius,
                const BoostBitset        & aLagrangeBasisActive,
                const uint               & aNumberOfActiveDesignLagrangeBasis,
                const Mat<real>          & aDimensions,
                const Mat<real>          & aDimensions_Offset,
                Mat< uint >              & aFoundBasis,
                Mat< real >              & aAlpha,
                Mat< real >              & aBeta
                );

                //Mat<real> const & aDimensionsOriginal,
                //Mat<real> const & aPointOfOrigin);
        /**
         * Provides a list of basis functions, which are in a specified square or cube on a certain level
         *
         * @param[in] aBasisID                       Basis function, which is used to search for neighbors
         * @param[in] aLevelToLook                   Level on which basis shall be searched
         * @param[in] aRadius                        Side length of bounding box
         * @param[in] aModelDim                      Dimension of Model
         * @param[in]  aNumberOfElementsPerDirection  As the name suggests
         *
         * @param[out] Mat             Output is a matrix containing basis IDs
         */
        void
        give_active_basis_on_level_within_bounding_box(
                const uint        & aBasisID,
                const uint        & aLevelToLook,
                const real        & aRadius,
                const Mat<real>   & aDimensions,
                const BoostBitset & aLagrangeBasisActive,
                const uint        & aModelDim,
                const uint        & aPolynomial,
                const Mat<uint>   & aNumberOfElementsPerDirection,
                Mat< uint >       & aBasisWithinBoundingBox,
                uint              & aBasisCounter );

        /**
         * subroutine for give_active_basis_on_level_within_bounding_box
         */
        void
        calculate_min_and_max_ijk(
                const uint        & aBasisID,
                const uint        & aLevelToLook,
                const real        & aRadius,
                const Mat<real>   & aDimensions,
                const uint        & aModelDim,
                const uint        & aPolynomial,
                const Mat<uint>   & aNumberOfElementsPerDirection,
                Mat < uint >      & aMinIJK,
                Mat < uint >      & aMaxIJK);

    };

}

#endif /* SRC_MESH_CL_LAGRANGE_FILTER_HPP_ */

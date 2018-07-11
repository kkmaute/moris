/*
 * cl_Hierarchical_Mesh_Filter.hpp
 *
 *  Created on: Jan 25, 2018
 *      Author: gleim
 */

#ifndef SRC_MESH_CL_HIERARCHICAL_MESH_FILTER_HPP_
#define SRC_MESH_CL_HIERARCHICAL_MESH_FILTER_HPP_

#include "algorithms.hpp"
#include "linalg.hpp"
#include "cl_Base_Mat.hpp" // LNA/src
#include "cl_Mat.hpp" // LNA/src
#include "cl_Cell.hpp" // CON/src
#include "fn_unique.hpp" // LNA/src
#include "fn_find_unique.hpp" // LNA/src
#include "fn_find.hpp" // LNA/src
#include "fn_sum.hpp" // LNA/src
#include "cl_BoostBitset.hpp" // CON/src

#include "cl_Hierarchical_Mesh_Basis.hpp"
#include "cl_Hierarchical_Mesh_Element.hpp"
#include "cl_Base_Mesh_Element.hpp"

namespace moris
{

    class Hierarchical_Mesh_Filter
    {
    protected:

    public:
        //Create Object of Basis
        Hierarchical_Mesh_Basis mBasis;
        //Create Object of Element
        Hierarchical_Mesh_Element mHMRElement;
        //Create Object of Element
        Base_Mesh_Element mBaseElement;
        /**
         * Hierarchical_Mesh constructor
         */
        Hierarchical_Mesh_Filter()
    {
    }

        /**
         * Hierarchical_Mesh destructor.
         */
        ~Hierarchical_Mesh_Filter() = default;

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
               BoostBitset & aDesignBSplineActive,
               Mat<real> const & aDimensions,
               Mat<real> const & aDimensions_Offset,
               Mat<real> const & aDimensionsOriginal,
               Mat<uint> & aIdFieldDesignNodalField,
               Mat<real> & aTMatrixDesignNodalField,
               bool const & aFilterLevelWeight,
               bool const & aPerformNormalization,
               Mat<real> const & aPointOfOrigin,
               Mat<uint> & IdFieldDesignNodalList);

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
        Mat<real>
        Filter_for_smoothing(
                uint const & aBasis,
                uint const & aDim,
                uint const & aPolynomialDesign,
                Mat<uint> const & aNumElements,
                uint const & aLevel,
                real const & aFilterRadius,
               BoostBitset & aElementActiveDesign,
               BoostBitset & aDesignBSplineActive,
               Mat<real> const & aDimensions,
               Mat<real> const & aDimensions_Offset);

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
        Mat<real>
        Filter_for_smoothing_new(
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
                Mat<real> const & aPointOfOrigin);

    };

}

#endif /* SRC_MESH_CL_HIERARCHICAL_MESH_FILTER_HPP_ */

/*
 * cl_Hierarchical_Mesh_Output.hpp
 *
 *  Created on: Jan 8, 2018
 *      Author: gleim
 */

#ifndef SRC_MESH_CL_HIERARCHICAL_MESH_OUTPUT_HPP_
#define SRC_MESH_CL_HIERARCHICAL_MESH_OUTPUT_HPP_

#include <iostream>
#include <sstream>
#include <fstream>

// MORIS libraries
#include "fn_sum.hpp" // LNA/src
#include "fn_save_vector_to_binary_file.hpp" // LNA/src
#include "cl_Map.hpp" // CON/src
#include "cl_Cell.hpp" // CON/src
#include "cl_Hierarchical_Mesh_MPI.hpp" // STK/src/Heirarchical
#include "cl_Hierarchical_Mesh_TMatrix.hpp" // STK/src/Heirarchical
//#include "cl_Hierarchical_Mesh_Basis.hpp" only needed for debugging // STK/src/Heirarchical

namespace moris
{

    class Hierarchical_Mesh_Output
    {
    protected:

    public:
// -----------------------------------------------------------------------------
        //Create Object of Basis
       // Hierarchical_Mesh_Basis mBasis;

        //Create Object of T-Matrix
        Hierarchical_Mesh_TMatrix mTMatrix;
// -----------------------------------------------------------------------------
        /**
         * Hierarchical_Mesh constructor
         */
        Hierarchical_Mesh_Output()
    {
    }

        /**
         * Hierarchical_Mesh destructor.
         */
        ~Hierarchical_Mesh_Output() = default;
// -----------------------------------------------------------------------------
        /**
         * Creates an output file for FemDoc with the T-Matrix and Id-Field
         *
         */
        void output_element_dof_connectivity_ascii(
                const uint              & aDim,
                const uint              & aPolynomial,
                const Mat<uint>         & aNumElements,
                const Mat<uint>         & aFetopo,
                const Mat<uint>         & aElementListOnProc,
                const BoostBitset       & aBasisActive,
                const bool              & aRefinement,
                const bool              & aTruncatedBsplines,
                Cell<Mat<real>>         & aTMatrixParentChild,
                Mat<uint>               & aBsplineBasisList,
                const map< uint, uint > & aBasisListMap );

// -----------------------------------------------------------------------------

        /**
         * Creates an output file for FemDoc with the T-Matrix and Id-Field
         *
         */
        void output_element_dof_connectivity_binary(
                const uint            & aDim,
                const uint            & aPolynomial,
                const Mat<uint>       & aNumElements,
                //const Mat<uint> & aNodalLocaltoGlobal,
                const Mat<uint>       & aFetopo,
                const Mat<uint>       & aElementListOnProc,
                const BoostBitset     & aBasisActive,
                const bool            & aRefinement,
                const bool            & aTruncatedBsplines,
                Cell<Mat<real>>       & aTMatrixParentChild,
                Mat<uint>             & aBsplineBasisList,
                const map< uint, uint > & aBasisListMap );

// -----------------------------------------------------------------------------

        /**
         * Creates an output file for FemDoc with the DOF variables and their interpolation per node
         *
         */
        void output_node_dof_connectivity_ascii(
                const Mat<uint>         & aLagrangeToBSplineMap,
                const Mat<uint>         & aIdFieldField,
                const Mat<real>         & aTMatrixField,
                const map< uint, uint > & aBasisListMap,
                const bool              & aRefinement);

// -----------------------------------------------------------------------------

        /**
         * Creates an output file for FemDoc with the DOF variables and their interpolation per node
         *
         */
        void output_node_dof_connectivity_binary(
                const Mat<uint>         & aLagrangeToBSplineMap,
                const Mat<uint>         & aIdFieldField,
                const Mat<real>         & aTMatrixField,
                const map< uint, uint > & aBasisListMap,
                const bool              & aRefinement);

// -----------------------------------------------------------------------------

        /**
         * Creates an output file for FemDoc with the design variables and their interpolation per node
         *
         */
        void output_design_variables_for_femdoc_ascii(
                const Mat<uint> & aLagrangeToBSplineMap,
                const Mat<uint> & aIdFieldFieldDesign,
                const Mat<real> & aTMatrixFieldDesign,
                const map<uint, uint> & aDesignBSplineListMap,
                const bool & aRefinement,
                const bool & aSymmetry
                );

// -----------------------------------------------------------------------------

        /**
         * Creates an output file for FemDoc with the design variables and their interpolation per node
         *
         */
        void output_design_variables_for_moris_ascii(
                const Mat<uint> & aLagrangeToBSplineMap,
                const Mat<uint> & aIdFieldFieldDesign,
                const Mat<real> & aTMatrixFieldDesign,
                const map<uint, uint> & aDesignBSplineListMap,
                const bool & aRefinement,
                const bool & aSymmetry
        );

// -----------------------------------------------------------------------------

        /**
         * Creates an output file for FemDoc with the design variables and their interpolation per node
         *
         */
        void output_design_variables_for_femdoc_binary(
                const Mat<uint> & aLagrangeToBSplineMap,
                const Mat<uint> & aIdFieldFieldDesign,
                const Mat<real> & aTMatrixFieldDesign,
                const map<uint, uint> & aDesignBSplineListMap,
                const bool & aRefinement,
                const bool & aSymmetry);

// -----------------------------------------------------------------------------

        /**
         * Creates an output file for FemDoc with the design variables and their interpolation per node
         *
         */
        void output_design_variables_for_moris_binary(
                const Mat<uint> & aLagrangeToBSplineMap,
                const Mat<uint> & aIdFieldFieldDesign,
                const Mat<real> & aTMatrixFieldDesign,
                const map<uint, uint> & aDesignBSplineListMap,
                const bool & aRefinement,
                const bool & aSymmetry);

// -----------------------------------------------------------------------------

        /**
         * Creates an output file for the continuation script to know in the next step the active design elements
         *
         */
        void output_active_design_elements(
                Mat<uint> aElementListActiveDesign,
                bool & aRefinement);

// -----------------------------------------------------------------------------

        /**
         * Creates an output file for the continuation script to know in the next step the active FEM elements
         *
         */
        void save_active_fem_elements(
                bool & aRefinement,
                Mat<uint> & aElementListOnProc);

// -----------------------------------------------------------------------------

        /**
         * Creates an output file for the continuation script to know in the next step the active design basis functions
         *
         */
        void save_active_design_basis(
                bool & aRefinement,
                Mat<uint> & aBasisActiveDesignList);
// -----------------------------------------------------------------------------

        /**
         * Creates an output file for the continuation script to know in the next step the active basis functions
         *
         */
        void save_active_basis(
                bool & aRefinement,
                Mat<uint> & aBasisActiveList);
// -----------------------------------------------------------------------------

        /**
         * Creates an output file for the continuation script to know in the next step the active basis functions
         *
         */
        void save_coordinate_list(
                const bool & aRefinement,
                const Mat<uint> & aNodalLocaltoGlobalExist);

// -----------------------------------------------------------------------------

        void save_field_to_file(
                const std::string      & aFilePath,
                const BoostBitset      & aBasisActive,
                const map< uint, real> & aFieldMap
                );

// -----------------------------------------------------------------------------

        static void
        make_design_variables_unique(
                Mat< uint > & aIdField,
                Mat< real > & aTMatrix );

// -----------------------------------------------------------------------------
    };
}

#endif /* SRC_MESH_CL_HIERARCHICAL_MESH_OUTPUT_HPP_ */

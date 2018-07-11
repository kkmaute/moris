/*
 * cl_Hierarchical_Mesh_Input.hpp
 *
 *  Created on: Jan 8, 2018
 *      Author: gleim
 */

#ifndef SRC_MESH_CL_HIERARCHICAL_MESH_INPUT_HPP_
#define SRC_MESH_CL_HIERARCHICAL_MESH_INPUT_HPP_

#include <iostream>
#include <sstream>
#include <fstream>
#include "cl_BoostBitset.hpp" // CON/src
#include "cl_Map.hpp" // CON/src
#include "linalg.hpp"
#include "cl_Base_Mat.hpp" // LNA/src
#include "cl_Mat.hpp" // LNA/src
#include "fn_load_vector_from_binary_file.hpp" // LNA/src
//#include "fn_save_vector_to_binary_file.hpp" // only needed for debugging // LNA/src
#include "cl_Base_Mesh_Element.hpp" // STK/src/Heirarchical
#include "cl_Hierarchical_Mesh_Basis.hpp" // STK/src/Heirarchical

namespace moris
{

    class Hierarchical_Mesh_Input
    {
    protected:

    public:
        //Create Object of Basis
        Hierarchical_Mesh_Basis mBasis;

        //Create Object of BaseElement
        Base_Mesh_Element mBaseElement;
        /**
         * Hierarchical_Mesh constructor
         */
        Hierarchical_Mesh_Input()
    {
    }

//-------------------------------------------------------------------------------

        /**
         * Hierarchical_Mesh destructor.
         */
        ~Hierarchical_Mesh_Input() = default;

//-------------------------------------------------------------------------------

        /**
         * Reads a list of active FEM elements from an older step (e.g. continuation script)
         *
         */
        void read_active_fem_elements(
                Mat<uint> & aElemLocaltoGlobalLastStepFEM,
                uint & aLevel,
                const uint & aDim,
                const Mat<uint> & aNumElements,
                BoostBitset & aElementActiveLastStep);

//-------------------------------------------------------------------------------

        /**
         * Reads a list of active FEM elements from an older step (e.g. continuation script)
         *
         */
        void read_active_design_elements(
                const uint      & aDim,
                const Mat<uint> & aNumElements,
                uint            & aLevel,
                Mat<uint>       & aElemLocaltoGlobalLastStepDesign,
                BoostBitset     & aElementActiveDesignLastStep);

//-------------------------------------------------------------------------------

        /**
         * Reads a list of active design basis functions from an older step (e.g. continuation script)
         *
         */
        void read_active_design_basis(
                const uint      & aDim,
                const uint      & aPolynomialDesign,
                const Mat<uint> & aNumElements,
                Mat<uint>       & aBasisActiveDesignListLastStep,
                BoostBitset     & aBasisActiveDesign);

//-------------------------------------------------------------------------------

        /**
         * Reads a list of active basis functions from an older step (e.g. continuation script)
         *
         */
        void read_active_basis(
                const uint      & aDim,
                const uint      & aPolynomialDesign,
                const Mat<uint> & aNumElements,
                Mat<uint>       & aBasisActiveDesignListLastStep,
                BoostBitset     & aBasisActive);

//-------------------------------------------------------------------------------

        /**
         * Reads the nodal ID design basis functions and T-Matrix from the designvariables.data from an older step ( e.g. continuation script)
         *
         */
        void read_design_variables_ascii(
                Mat<uint> & aIdFieldDesignNodalField,
                Mat<real> & aTMatrixDesignNodalField);

//-------------------------------------------------------------------------------

        /**
         * Reads the nodal ID design basis functions and T-Matrix from the designvariables.data from an older step ( e.g. continuation script)
         *
         */
        void read_design_variables_binary(
                Mat<uint> & aIdFieldDesignNodalField,
                Mat<real> & aTMatrixDesignNodalField);

//-------------------------------------------------------------------------------

        /**
         * Reads a list of coordinate numbers of the last step
         *
         */
        void read_coordinate_list(
                 Mat<uint> & aCoordinateList);

//-------------------------------------------------------------------------------

        /**
         * Reads a field of signed distances from a file
         *
         */
        void read_sdf_file(
                 std::string const & aSDFFileName,
                 Mat<real> & aSDFNodalField,
                 Mat<uint> & aCoordinateList);

//-------------------------------------------------------------------------------

        /**
         * Reads data from the AbsDesVariables file and puts user specific calculations on the different fields
         *
         */
        void read_absdesvariables_file(
                const uint                      & aDim,
                const uint                      & aPolynomialDesign,
                const Mat<uint>                 & aNumElements,
                const Mat<uint>                 & aNumBasis,
                const Mat<uint>                 & aBasisActiveDesignListLastStep,
                map<uint, real>                 & aNodalADV,
                map<uint, real>                 & aNodalPDV,
                map<uint, real>                 & aNodalLvlsetField,
                const Cell< map< uint, real > > & aNodalSDFs,
                BoostBitset                     & aNodalFieldExists,
                Mat<uint>                       & tIdFieldDesignNodalField,
                Mat<real>                       & tTMatrixDesignNodalField,
                real                            & tSimpExp,
                real                            & sElementEdgeWidth,
                real                            & sInitialElementEdgeWidth,
                real                            & sLSscale,
                real                            & sLSthresh,
                real                            & sDensityShift,
                real                            & sPrjbeta,
                const Mat< real >               & aSDFAlpha,
                const Mat< real >               & aSDFBeta,
                const bool                      & aExpectBinaryInput,
                const bool                      & aUseSymmetry,
                const uint                      & aSymmetryPlane,
                const uint                      & aSymmetryIndex );

//-------------------------------------------------------------------------------

        /**
          * needed to create nodal field if AbsDesVariables.dat is symmetric
          */
        void expand_symmetric_field(
                        Mat< real >       & aField,
                        const Mat< uint > & aBasisActiveDesignListLastStep,
                        const uint        & aModelDim,
                        const uint        & aPolynomial,
                        const uint        & aSymmetryPlane,
                        const uint        & aSymmetryIndex,
                        const Mat<uint>   & aNumberOfElementsPerDirection,
                        const Mat<uint>   & aNumberOfBasisPerDirection );

//-------------------------------------------------------------------------------

        /**
         * Reads a nodal field from a bianry file and creates a map
         */
        void
        read_field_from_file(
                const std::string & aFilePath,
                const BoostBitset & aBasisActive,
                const bool        & aUseSymmetry,
                const Mat< uint > & aBasisActiveDesignListLastStep,
                const uint        & aModelDim,
                const uint        & aPolynomial,
                const uint        & aSymmetryPlane,
                const uint        & aSymmetryIndex,
                const Mat<uint>   & aNumberOfElementsPerDirection,
                const Mat<uint>   & aNumberOfBasisPerDirection,
                map< uint, real>  & aFieldMap);



    };
}
#endif /* SRC_MESH_CL_HIERARCHICAL_MESH_INPUT_HPP_ */

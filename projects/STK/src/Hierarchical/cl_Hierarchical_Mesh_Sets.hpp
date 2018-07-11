/*
 * cl_Hierarchical_Mesh_Sets.hpp
 *
 *  Created on: Jan 5, 2018
 *      Author: gleim
 */

#ifndef SRC_MESH_CL_HIERARCHICAL_MESH_SETS_HPP_
#define SRC_MESH_CL_HIERARCHICAL_MESH_SETS_HPP_

#include "linalg.hpp"
#include "cl_Base_Mat.hpp" // LNA/src
#include "cl_Mat.hpp" // LNA/src
#include "fn_find.hpp" // LNA/src
#include "fn_sum.hpp" // LNA/src
#include "cl_Cell.hpp" // CON/src
#include "cl_BoostBitset.hpp" // CON/src
#include "cl_Hierarchical_Mesh_Element.hpp"
#include "cl_Hierarchical_Mesh_Basis.hpp"
#include "cl_Hierarchical_Mesh_MPI.hpp"
#include "cl_Base_Mesh_Element.hpp"

namespace moris
{

    class Hierarchical_Mesh_Sets
    {
    protected:

    public:
        //Create Object of Basis
        Hierarchical_Mesh_Basis mBasis;
        //Create Object of Element
        Hierarchical_Mesh_Element mHMRElement;
        //Create Object of MPI
        Hierarchical_Mesh_MPI mMPI;
        //Create Object of base Element
        Base_Mesh_Element mBaseElement;
        /**
         * Hierarchical_Mesh constructor
         */
        Hierarchical_Mesh_Sets()
    {
    }

        /**
         * Hierarchical_Mesh destructor.
         */
        ~Hierarchical_Mesh_Sets() = default;

        /**
         * Creates a distributed nodeset ( a nodeset can be defined and the respective processor grabs the nodes which it owns)
         *
         * @param[in] aModelDim                     Dimensions of the problem
         * @param[in] aPolynomial                   Polynomial degree of the problem
         * @param[in] aNumberOfElementsPerDirection   Number of elements in each direction
         * @param[in] aInNodeSet                    General nodeset which is user defined
         * @param[in] aElementListOnProc            List of elements, which are active on this proc
         *
         * @param[out] aOutNodeSet   Distributed nodeset for each proc
         *
         */
        Mat<uint>
        set_nodeset(
                uint const & aModelDim,
                uint const & aPolynomial,
                Mat<uint> const & aNumberOfElementsPerDirection,
                Mat<uint> const & aInNodeSet,
                Mat<uint> const & aElementListOnProc);

        /**
         * Updates the NodeSet for the active elements/ basis functions (Only possible if aPolynomial = 1, otherwise the B-spline coefficients never overlap ! )
         *
         * @param[in] aModelDim                     Dimensions of the problem
         * @param[in] aPolynomial                   Polynomial degree of the problem
         * @param[in] aNumberOfElementsPerDirection   Number of elements in each direction
         * @param[in] aLevel                        Current level (level of finest elements)
         * @param[in] aNumSets                      Number of node sets
         * @param[in] aInNodeSet                    General nodeset which is user defined
         * @param[in] aBasisActive                  Bitset with active basis functions
         *
         * @param[out] aNodeSet                     Nodest is updated from the initial to the current level
         *
         */
        Cell<Mat<uint>>
        update_nodeset(
                uint const & aModelDim,
                uint const & aPolynomial,
                Mat<uint> const & aNumberOfElementsPerDirection,
                uint const & aLevel,
                Mat<uint> const & aNumSets,
                Mat<uint> & aNodeSet,
                BoostBitset const & aBasisActive);

        /**
         * Creates a distributed side set. Each processor takes from the general side set the part he owns
         *
         * @param[in] aModelDim                     Dimensions of the problem
         * @param[in] aPolynomial                   Polynomial degree of the problem
         * @param[in] aNumberOfElementsPerDirection   Number of elements in each direction
         * @param[in] aDecomp                       Decomposition of the mesh (Processor needs to know his domain)
         * @param[in] aSideSet                      General sideset for the whole domain (Each processor grabs his necessary part)
         *
         * @param[out] aOutSideSet                  Distributed side set for each proc
         *
         */
        Mat<uint>
        set_sideset(
                uint const & aModelDim,
                uint const & aPolynomial,
                Mat<uint> const & aNumberOfElementsPerDirection,
                Mat<uint> const & aDecomp,
                Mat<uint> const & aSideSet);

        /**
         * Updates the SideSet for the active elements/ basis functions
         *
         * @param[in] aModelDim                     Dimensions of the problem
         * @param[in] aNumberOfElementsPerDirection   Number of elements in each direction
         * @param[in] aLevel                        Finest level of refinement
         * @param[in] aNumSets                      Number of side sets  (Mat<uint> includes number of node, side, block sets)
         * @param[in] aSideSet                      Distributed side set for each proc
         * @param[in] aElementActive                Active elements in a bitset to know which side sets are needed
         *
         * @param[out] aOutSideSet                  Cell with all side sets (in each cell is a side set)
         *
         */
        Cell<Mat<uint>>
        update_sideset(
                uint const & aModelDim,
                Mat<uint> const & aNumberOfElementsPerDirection,
                uint const & aLevel,
                Mat<uint> const & aNumSets,
                Mat<uint> const & aSideSet,
                BoostBitset const & aElementActive);

        /**
         * Updates the BlockSet for the active elements
         *
         * @param[in] aModelDim                     Dimensions of the problem
         * @param[in] aNumberOfElementsPerDirection   Number of elements in each direction
         * @param[in] aNumSets                      Number of block sets  (Mat<uint> includes number of node, side, block sets)
         * @param[in] aBlockSet                     Distributed block set for each proc
         * @param[in] aElementListOnProc            List of active elements
         *
         * @param[out] aOutBlockSet                 Cell with all block sets (in each cell is a side set)
         *
         */
        Cell<Mat<uint>>
        update_blockset(
                uint const & aModelDim,
                Mat<uint> const & aNumberOfElementsPerDirection,
                Mat<uint> const & aNumSets,
                Mat<uint> const & aBlockset,
                Mat<uint> const & aElementListOnProc);

    };
}

#endif /* SRC_MESH_CL_HIERARCHICAL_MESH_SETS_HPP_ */

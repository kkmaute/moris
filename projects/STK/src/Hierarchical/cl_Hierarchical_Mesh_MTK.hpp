/*
 * cl_Hierarchical_Mesh_MTK.hpp
 *
 *  Created on: Jan 4, 2018
 *      Author: gleim
 */

#ifndef SRC_MESH_CL_HIERARCHICAL_MESH_MTK_HPP_
#define SRC_MESH_CL_HIERARCHICAL_MESH_MTK_HPP_

#include "cl_Mesh.hpp" // MTK/src
#include "cl_Cell.hpp" // CON/src
#include "linalg.hpp"
#include "cl_Base_Mesh_Element.hpp"

namespace moris
{

    class Hierarchical_Mesh_MTK
    {
    protected:

    public:
        //Create Object of BaseElement
        Base_Mesh_Element mBaseElement;
        /**
         * Hierarchical_Mesh constructor
         */
        Hierarchical_Mesh_MTK()
    {
    }

        /**
         * Hierarchical_Mesh destructor.
         */
        ~Hierarchical_Mesh_MTK() = default;

        /**
         * Creates the MTK database to output a STK file
         *
         * @param[in] OutputFileName                Name of the output file
         * @param[in] aModelDim                     Number of dimensions.
         * @param[in] aNumberOfElementsPerDirection   Number of elements in each direction =[Number of elements in x,y,z direction].
         * @param[in] aFetopo                       Topology of the elements: rows = elements, cols = connected nodes
         * @param[in] aFeCoord                      Coordinates
         * @param[in] aBasisProcOwner               A list of ownership of basis/nodes
         * @param[in] aElementListOnProc            A list of active elements
         * @param[in] aLagrangeToBSplineMap         Map of nodal information (includes hanging nodes with different numbers)
         * @param[in] aLagrangeListOnProc           Map of nodal information with original basis id's
         * @param[in] aSideSet                      Side sets of the domain
         *
         */
        void
        create_MTK_file(
                std::string const & OutputFileName,
                uint & aModelDim,
                Mat<uint> & aNumberOfElementsPerDirection,
                Mat<uint> & aFetopo,
                Mat<real> & aFeCoord,
                Mat<uint> & aBasisProcOwner,
                Mat<uint> & aElementListOnProc,
                Mat<uint> & aLagrangeToBSplineMap,
                Mat<uint> & aLagrangeListOnProc,
                Cell<Mat<uint>> & aSideSet);

        /**
         * Creates the MTK database to output a STK file with field data
         *
         * @param[in] OutputFileName                Name of the output file
         * @param[in] aModelDim                     Number of dimensions.
         * @param[in] aNumberOfElementsPerDirection   Number of elements in each direction =[Number of elements in x,y,z direction].
         * @param[in] aFetopo                       Topology of the elements: rows = elements, cols = connected nodes
         * @param[in] aFeCoord                      Coordinates
         * @param[in] aBasisProcOwner               A list of ownership of basis/nodes
         * @param[in] aElementListOnProc            A list of active elements
         * @param[in] aLagrangeToBSplineMap           Map of nodal information (includes hanging nodes with different numbers)
         * @param[in] aLagrangeListOnProc       Map of nodal information with original basis id's
         * @param[in] aSideSet                      Side sets of the domain
         * @param[in] aFieldData                    Field data for MTK for outputing elemental or nodal data
         * @param[in] aFieldName                    Name of the fiels
         * @param[in] aFieldRank                    Entity rank of the fields
         *
         */
        void create_MTK_file(
                std::string const & OutputFileName,
                uint & aModelDim,
                Mat<uint> & aNumberOfElementsPerDirection,
                Mat<uint> & aFetopo,
                Mat<real> & aFeCoord,
                Mat<uint> & aBasisProcOwner,
                Mat<uint> & aElementListOnProc,
                Mat<uint> & aLagrangeToBSplineMap,
                Mat<uint> & aLagrangeListOnProc,
                Cell<Mat<uint>> & aSideSet,
                Cell<Mat<real>> & aFieldData,
                Cell<std::string> & aFieldName,
                Cell<enum EntityRank> & aFieldRank);

    };
}

#endif /* SRC_MESH_CL_HIERARCHICAL_MESH_MTK_HPP_ */

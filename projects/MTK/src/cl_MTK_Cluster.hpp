/*
 * cl_MTK_Cluster.hpp
 *
 *  Created on: Aug 17, 2018
 *      Author: messe
 */

#ifndef PROJECTS_MTK_SRC_CL_MTK_CLUSTER_HPP_
#define PROJECTS_MTK_SRC_CL_MTK_CLUSTER_HPP_

#include "typedefs.hpp" //MRS/COR/src
#include "cl_Cell.hpp" //MRS/CON/src
#include "cl_Mat.hpp" //LNA/src
#include "cl_MTK_Vertex.hpp" //MTK/src
#include "cl_MTK_Enums.hpp" //MTK/src
#include "cl_MTK_Cell.hpp" //MTK/src

//------------------------------------------------------------------------------
namespace moris
{
    namespace mtk
    {
//------------------------------------------------------------------------------

    class Cluster
    {
//------------------------------------------------------------------------------
    public:
//------------------------------------------------------------------------------

        /**
         * trivial constructor
         */
        Cluster(){};

//------------------------------------------------------------------------------

        /**
         * Destructor. Must be virtual.
         */
        virtual
        ~Cluster(){};

//------------------------------------------------------------------------------

        /**
         * returns the number of cells in this cluster
         *
         * @return uint number of cells
         */
        virtual uint
        get_number_of_cells() const = 0;

//------------------------------------------------------------------------------

        /**
         * returns the cell pointers of this cluster
         *
         * @return uint number of cells
         */
        virtual moris::Cell< mtk::Cell* >
        get_cell_pointers() = 0;

//------------------------------------------------------------------------------

    };
//------------------------------------------------------------------------------
    } /* namespace mtk */
} /* namespace moris */
//------------------------------------------------------------------------------


#endif /* PROJECTS_MTK_SRC_CL_MTK_CLUSTER_HPP_ */

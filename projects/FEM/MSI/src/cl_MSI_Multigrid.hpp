/*
 * cl_MSI_Multigrid.hpp
 *
 *  Created on: Jul 14, 2018
 *      Author: schmidt
 */
#ifndef SRC_FEM_CL_MSI_MULTIGRID_HPP_
#define SRC_FEM_CL_MSI_MULTIGRID_HPP_

#include "cl_Cell.hpp"
#include "cl_Communication_Tools.hpp"
#include "cl_Communication_Manager.hpp"
#include "cl_MSI_Model_Solver_Interface.hpp"
#include "cl_Map.hpp"
#include "fn_sum.hpp"

namespace moris
{
    namespace MSI
    {
        class Multigrid
        {
        private:

        public:
            Multigrid(){};

            ~Multigrid(){};

        };
    } /* namespace MSI */
} /* namespace moris */

#endif /* SRC_FEM_CL_MSI_MULTIGRID_HPP_ */

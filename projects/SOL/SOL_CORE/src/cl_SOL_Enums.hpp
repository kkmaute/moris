/*
 * cl_SOL_Enums.hpp
 *
 *  Created on: Jul 2, 2018
 *      Author: schmidt
 */
#ifndef SRC_DISTLINALG_CL_DISTLINALG_ENUMS_HPP_
#define SRC_DISTLINALG_CL_DISTLINALG_ENUMS_HPP_

namespace moris
{
    namespace sol
    {
        enum class SolverType
        {
            TRILINOSTEST,     //< Wrapper around Aztec Solver
            AZTEC_IMPL,            //< Wrapper around Aztec Solver
            AMESOS_IMPL,           //< Wrapper around Amesos Solver
            AMESOS2_IMPL,          //< Wrapper around Amesos2 Solver
            BELOS_IMPL,          //< Wrapper around Belos Solver
            PETSC             //< Wrapper around Petsc Solver
        };

        enum class MapType
        {
            Epetra,   // Indicates the Vector/Matrix/Map type
            Petsc     // Indicates the Vector/Matrix/Map type
        };
    }
}

#endif /* SRC_DISTLINALG_CL_DISTLINALG_ENUMS_HPP_ */

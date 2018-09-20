/*
 * cl_Model_Solver_Interface_Solver.hpp
 *
 *  Created on: Mar 28, 2018
 *      Author: schmidt
 */
#ifndef SRC_DISTLINALG_CL_MODEL_SOLVER_INTERFACE_HPP_
#define SRC_DISTLINALG_CL_MODEL_SOLVER_INTERFACE_HPP_

// MORIS header files.
#ifdef MORIS_HAVE_PARALLEL
#include <mpi.h>
#endif

#include <iostream>
#include <chrono>
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

#include "cl_Matrix_Vector_Factory.hpp" // DLA/src
#include "cl_Linear_Solver.hpp"         // DLA/src
#include "cl_Sparse_Matrix.hpp"         // DLA/src
#include "cl_Vector.hpp"         // DLA/src
#include "cl_Solver_Input.hpp"

typedef std::chrono::high_resolution_clock Clock;

namespace moris
{
//class Solver_Input;
class Model_Solver_Interface
{
private:
protected:

public:

    Model_Solver_Interface( moris::Linear_Solver * aLin,
                            moris::Solver_Input  * aInput,
                            moris::Sparse_Matrix * aMat,
                            moris::Dist_Vector   * aVectorRHS );

//---------------------------------------------------------------------------------------------------------
    ~Model_Solver_Interface();
};
}

#endif /* SRC_DISTLINALG_CL_MODEL_SOLVER_INTERFACE_HPP_ */

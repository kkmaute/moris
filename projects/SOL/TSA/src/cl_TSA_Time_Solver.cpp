/*
 * cl_TSA_Time_Solver.cpp
 *
 *  Created on: Feb 02, 2019
 *      Author: schmidt
 */
#include <ctime>

#include "cl_Vector.hpp"

#include "cl_Communication_Tools.hpp"
#include "cl_TSA_Time_Solver.hpp"

using namespace moris;
using namespace tsa;

//-------------------------------------------------------------------------------

Time_Solver::Time_Solver()
{
    std::cout<<"Time Solver"<<std::endl;
}

//-------------------------------------------------------------------------------

void Time_Solver::get_full_solution( moris::Matrix< DDRMat > & LHSValues )
{
    mFullVector->extract_copy( LHSValues );
}


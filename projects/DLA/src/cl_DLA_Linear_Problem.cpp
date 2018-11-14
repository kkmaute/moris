/*
 * cl_DLA_Linear_Problem.cpp
 *
 *  Created on: Dec 6, 2017
 *      Author: schmidt
 */
#include "cl_DLA_Linear_Problem.hpp"
#include "cl_Vector.hpp"

namespace moris
{
namespace dla
{
Dist_Vector * Linear_Problem::get_full_solver_LHS()
{
    mFullVectorLHS->vec_put_scalar( 0.0 );

    mFullVectorLHS->import_local_to_global( *mFreeVectorLHS );

    return mFullVectorLHS;
}

void Linear_Problem::set_free_solver_LHS( Dist_Vector * aFullSolVector)
{
    mFreeVectorLHS->import_local_to_global( *aFullSolVector );
}

}
}


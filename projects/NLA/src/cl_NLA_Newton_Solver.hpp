/*
 * cl_NLA_Newton_Solver.hpp
 *
 *  Created on: Sep 21, 2018
 *      Author: schmidt
 */
#ifndef SRC_FEM_CL_NEWTON_SOLVER_HPP_
#define SRC_FEM_CL_NEWTON_SOLVER_HPP_

#include "typedefs.hpp"

#include "typedefs.hpp"

#include "cl_Matrix_Vector_Factory.hpp"

namespace moris
{
class Dist_Vector;
class Linear_Solver;
namespace NLA
{
    class Newton_Solver
    {
    private:
        moris::uint mA;

        Dist_Vector   * mVectorFull;
        Dist_Vector   * mVectorMaster;
        Dist_Vector   * mVectorMasterSol;
        Dist_Vector   * mVectorMasterPrev;

        Map_Class     * mMap;

        std::shared_ptr< Linear_Solver > mLinearSolver;

    public:
        Newton_Solver()
        {
        };

        Newton_Solver( std::shared_ptr< Linear_Solver > aLinearSolver );

        ~Newton_Solver()
        {};

        void devide( moris::uint aA )
        {
            mA = mA/aA;
        };

    };

}
}



#endif /* SRC_FEM_CL_NEWTON_SOLVER_HPP_ */

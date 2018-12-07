/*
 * cl_DLA_Geometric_Multigrid.hpp
 *
 *  Created on: Apr 6, 2018
 *      Author: schmidt
 */

#ifndef SRC_DISTLINALG_CL_DLA_GEOMETRIC_MULTIGRID_HPP_
#define SRC_DISTLINALG_CL_DLA_GEOMETRIC_MULTIGRID_HPP_

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

namespace moris
{
    class Dist_Vector;
    class Sparse_Matrix;
    class Solver_Interface;

    namespace mtk
    {
        class Mesh;
    }
    namespace dla
    {
    class Geometric_Multigrid
    {
    private:
        Solver_Interface * mSolverInterface;

        mtk::Mesh * mMesh;

        moris::Cell< Matrix< DDUMat > > mListAdofExtIndMap;

    public:
        Geometric_Multigrid( Solver_Interface * aSolverInterface );

        /** Destructor */
        ~Geometric_Multigrid(){};





        virtual mtk::Mesh * get_mesh_pointer_for_multigrid()
        {
            MORIS_ERROR(false, "Solver_Interface::get_mesh_pointer_for_multigrid, Only works with MSI and multigrid");
            return nullptr;
        };


    };
}
}


#endif /* SRC_DISTLINALG_CL_DLA_GEOMETRIC_MULTIGRID_HPP_ */

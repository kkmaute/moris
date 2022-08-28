/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_DLA_Geometric_Multigrid.hpp
 *
 */

#ifndef SRC_DISTLINALG_CL_DLA_GEOMETRIC_MULTIGRID_HPP_
#define SRC_DISTLINALG_CL_DLA_GEOMETRIC_MULTIGRID_HPP_

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

namespace moris
{
    namespace sol
    {
        class Dist_Vector;
        class Dist_Matrix;
    }
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
        //! List containing the prolongation operators
        moris::Cell< sol::Dist_Matrix * > mProlongationList;

        //! List of external indices for each level
        moris::Cell< Matrix< DDUMat > > mListAdofExtIndMap;

        //! List of type/time identifiers for each level
        moris::Cell< Matrix< DDSMat > > mListAdofTypeTimeIdentifier;

        //! Map which maps external indices to internal MSI indices. List 1 = Level; List 2 = type/time;
        moris::Cell< moris::Cell< Matrix< DDSMat > > > mMultigridMap;

        // Pointer to model solver interface
        Solver_Interface * mSolverInterface;

        // Mesh pointer
        mtk::Mesh * mMesh;

    public:
        /**
         * @brief Constructor. Build the list with prolongation operators
         *
         * @param[in] aSolverInterface    Pointer to solverInterface
         *
         */
        Geometric_Multigrid( Solver_Interface * aSolverInterface );

        /** Destructor */
        ~Geometric_Multigrid();

        /**
         * @brief Returns list with operators
         *
         * @param[out] mProlongationList    List with prolongation operators
         *
         */
        moris::Cell< sol::Dist_Matrix * > & get_prolongation_list()
        {
            return mProlongationList;
        };
    };
}
}

#endif /* SRC_DISTLINALG_CL_DLA_GEOMETRIC_MULTIGRID_HPP_ */


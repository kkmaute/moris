/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_EqnManager.hpp
 *
 */

#ifndef MORIS_MODEL_CL_EQNMANAGER_HPP_
#define MORIS_MODEL_CL_EQNMANAGERL_HPP_

// MORIS project header files.
//#include "cl_Mesh_Mtk.hpp"
#include "core.hpp"
#include "cl_Cell.hpp" // CON/src
#include "ios.hpp"
#include "cl_Model_Enums.hpp" // MOD/src

// MORIS library header files.
#include "algorithms.hpp"
#include "core.hpp"
#include "linalg.hpp"
#include "cl_Mat.hpp"

namespace moris
{
    // Argument list for build_global_system
    struct ArgListEqnMgr_build_global_system
    {
        moris::Mat< moris::real > GlbSolVecU;
        moris::Mat< moris::real > GlbSolVecV;
        moris::Mat< moris::real > GlbSolVecA;
        moris::Mat< moris::real > ActiveDoFs;
        moris::real*              TimeIntegrator;
        moris::real*              Solver;
    };

    // Argument list for get_eqn_obj
    struct ArgListEqnMgr_get_eqn_obj
    {

        enum EquationObjectType   EqnObjType;
        moris::Mat< moris::real > GlbSolVecU;
        moris::Mat< moris::real > ResVec;
        moris::Mat< moris::real > MMat;
        moris::Mat< moris::real > CMat;
        moris::Mat< moris::real > KMat;
        moris::uint               EqnObjId;
    };

    class EqnManager
    {
    protected:

        moris::uint mNumberOfEqnObjects = 0; // Default value of number of equation objects member variable

    public:

        moris::ArgListEqnMgr_build_global_system mArgListEqnMgr_build_global_system; // Build_global_system specific argument list
        moris::ArgListEqnMgr_get_eqn_obj         mArgListEqnMgr_get_eqn_obj;         // Get_eqn_obj specific argument list

        /**
         * EqnManager constructor
         */
        EqnManager()
        {

        }

        /**
         * Eqn_Manager destructor.
         */
        ~EqnManager() = default;

        /**
         * Creates based on the global solution vector a residual vector and a jacobian matrix and solves this system of equations.
         *
         * @param[in] aL       .... Argument list for build_global_system.
         * @param[out] aL      .... Argument list for build_global_system.
         *
         */
        void build_global_system(moris::ArgListEqnMgr_build_global_system & aL);

        /**
         * Calculates the residual vector, the sparse jacobian matrix, the row and coloumn vector of the sparse jacobian matrix.
         *
         * @param[in] aL       .... Argument list for get_eqn_obj.
         * @param[out] aL      .... Argument list for get_eqn_obj.
         *
         */
        void get_eqn_obj(moris::ArgListEqnMgr_get_eqn_obj & aL);

        /**
         * Checks the residual vector, jacobian matrix, the col and row vectors for the sparsity pattern
         *
         * @param[in] aResVec     .... Residual vector.
         * @param[in] aJacMat     .... Jacobian matrix in a vector (sparsity pattern).
         * @param[in] aJacMat_col .... Coloumn vector of the sparse jacboian matrix
         * @param[in] aJacMat_row .... Row vector of the sparse jacboian matrix
         *
         *  @note The residual vector, the jacobian matrix, the row and column vector of the sparse jacobian matrix are checked for finite entries.\n
         *        In addition, the vectors are checked for positive entries.
         *
         */
        void check_eqn_obj(moris::Mat<moris::real> &aResVec,
                moris::Mat<moris::real> &aJacMat,
                moris::Mat<moris::uint> &aJacMat_col,
                moris::Mat<moris::uint> &aJacMat_row,
                moris::uint &aEqnObjId);

        /**
         * Set number of equation objects as a member variable of Eqn_Manager class.
         *
         * @param[in] aNumberOfEqnObjects  .... Number of equation objects.
         *
         */
        void set_NumberOfEqnObjects(moris::real aNumberOfEqnObjects)
        {
            mNumberOfEqnObjects = aNumberOfEqnObjects;
        }
    };

}   // namespace moris

#endif /* MORIS_MODEL_CL_EQNMANAGER_HPP_ */


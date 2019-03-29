/*
 * cl_NLA_Nonlinear_Problem.hpp
 *
 *  Created on: Nov 18, 2018
 *      Author: schmidt
 */
#ifndef MORIS_DISTLINALG_CL_NLA_NONLINEAR_PROBLEM_HPP_
#define MORIS_DISTLINALG_CL_NLA_NONLINEAR_PROBLEM_HPP_

// MORIS header files.
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

#include "cl_NLA_Nonlinear_Solver_Enums.hpp"
#include "cl_DLA_Linear_Problem.hpp"

#include "cl_Param_List.hpp"

namespace moris
{
class Map_Class;
class Dist_Vector;
class Solver_Interface;
namespace dla
{
    class Linear_Problem;
}
namespace NLA
{
    class SOL_Warehouse;
    class Nonlinear_Problem
    {
    private:

        void  delete_pointers();

        //--------------------Arc Length-------------------
        Sparse_Matrix * mFullForKTilde  = nullptr;

        Dist_Vector * mFullForKDiag     = nullptr;
        Dist_Vector * mFullForK0Diag    = nullptr;

        Dist_Vector * mDTildeVec        = nullptr;
        Dist_Vector * mDTilde0Vec       = nullptr;

        Dist_Vector * mFext             = nullptr;
        //-------------------------------------------------

    protected:
        Dist_Vector * mFullVector = nullptr;

        Map_Class   * mMap = nullptr;
        Map_Class   * mMapFull = nullptr;               //FIXME replace with marketplace

        dla::Linear_Problem * mLinearProblem = nullptr;

        bool mIsMasterSystem = false;

        bool mBuildLinerSystemFlag = true;

        //! Map type. for special Petsc functionalities
        enum MapType mMapType = MapType::Epetra;

        //! Nonlinear solver manager index. only for output purposes
        moris::sint mNonlinearSolverManagerIndex = -1;


    public:

        //--------------------------------------------------------------------------------------------------
        /**
         * @brief Constructor. Creates nonlinear system
         *
         * @param[in] aSolverInterface             Pointer to the solver interface
         * @param[in] aNonlinearSolverManagerIndex Nonlinera solver manager index. Default = 0
         * @param[in] aBuildLinerSystemFlag        Flag if linear system shall be build or not. Default = true
         * @param[in] aMapType                     Map type. Epetra or Petsc. Default MapType::Epetra
         */
        Nonlinear_Problem(       Solver_Interface * aSolverInterface,
                           const moris::sint        aNonlinearSolverManagerIndex = 0,
                           const bool               aBuildLinerSystemFlag = true,
                           const enum MapType       aMapType = MapType::Epetra );

        //--------------------------------------------------------------------------------------------------
        /**
         * @brief Constructor. Creates nonlinear system
         *
         * @param[in] aNonlinDatabase             Pointer to database
         * @param[in] aSolverInterface             Pointer to the solver interface
         * @param[in] aNonlinearSolverManagerIndex Nonlinera solver manager index. Default = 0
         * @param[in] aBuildLinerSystemFlag        Flag if linear system shall be build or not. Default = true
         * @param[in] aMapType                     Map type. Epetra or Petsc. Default MapType::Epetra
         */
        Nonlinear_Problem(       SOL_Warehouse    * aNonlinDatabase,
                                 Solver_Interface * aSolverInterface,
                                 Dist_Vector      * aFullVector,
                           const moris::sint        aNonlinearSolverManagerIndex = 0,
                           const bool               aBuildLinerSystemFlag = true,
                           const enum MapType       aMapType = MapType::Epetra);

        //--------------------------------------------------------------------------------------------------
        ~Nonlinear_Problem();

        //--------------------------------------------------------------------------------------------------
        void set_interface( Solver_Interface * aSolverInterface );

        //--------------------------------------------------------------------------------------------------
        void build_linearized_problem( const bool & aRebuildJacobian,
                                             sint   aNonLinearIt );

        //--------------------------------------------------------------------------------------------------
        void build_linearized_problem( const bool & aRebuildJacobian,
                                       const sint   aNonLinearIt,
                                       const sint   aRestart );

        //--------------------------------------------------------------------------------------------------
        void print_sol_vec( const sint aNonLinearIt );

        //--------------------------------------------------------------------------------------------------
        void restart_from_sol_vec( const sint aNonLinearIt );

        //--------------------------------------------------------------------------------------------------
        dla::Linear_Problem * get_linearized_problem(){ return mLinearProblem; };


        //--------------------------------------------------------------------------------------------------
        Dist_Vector * get_full_vector();

        //--------------------------------------------------------------------------------------------------
        void extract_my_values( const moris::uint             & aNumIndices,
                                const moris::Matrix< DDSMat > & aGlobalBlockRows,
                                const moris::uint             & aBlockRowOffsets,
                                      moris::Matrix< DDRMat > & LHSValues );

        //--------------------------------------------------------------------------------------------------
        //--------------------------------Arc-Length 'get' functions----------------------------------------
        Dist_Vector * get_k_diag()
        {
            return mFullForKDiag;
        }
        //--------------------------------------------------------------------------------------------------
        Dist_Vector * get_k0_diag()
        {
            return mFullForK0Diag;
        }
        //--------------------------------------------------------------------------------------------------
        Sparse_Matrix * get_full_for_consistent_tangent()
        {
            return mFullForKTilde;
        }
        //--------------------------------------------------------------------------------------------------
        Dist_Vector * get_d_tilde()
        {
            return mDTildeVec;
        }
        //--------------------------------------------------------------------------------------------------
        Dist_Vector * get_d_tilde0()
        {
            return mDTilde0Vec;
        }
        //--------------------------------------------------------------------------------------------------
        Dist_Vector * get_f_ext()
        {
            return mFext;
        }
    };
}
}
#endif /* MORIS_DISTLINALG_CL_NLA_NONLINEAR_PROBLEM_HPP_ */

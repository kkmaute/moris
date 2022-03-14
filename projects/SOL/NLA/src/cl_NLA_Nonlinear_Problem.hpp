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
#include "cl_SOL_Enums.hpp"

#include "cl_Param_List.hpp"

namespace moris
{
    namespace sol
    {
        class Dist_Map;
        class Dist_Vector;
    }    // namespace sol
    class Solver_Interface;
    namespace dla
    {
        class Linear_Problem;
    }
    namespace sol
    {
        class SOL_Warehouse;
    }
    namespace NLA
    {
        class Nonlinear_Solver;
        class Nonlinear_Problem
        {
          private:
            void delete_pointers();

          protected:
            Solver_Interface* mSolverInterface;

            sol::Dist_Vector* mFullVector      = nullptr;
            sol::Dist_Vector* mDummyFullVector = nullptr;    // FIXME Delete

            sol::Dist_Map* mMap     = nullptr;
            sol::Dist_Map* mMapFull = nullptr;    // FIXME replace with marketplace

            dla::Linear_Problem* mLinearProblem = nullptr;

            //! Pointer to solver database
            sol::SOL_Warehouse* mSolverWarehouse = nullptr;

            //! Pointer to my nonlinear solver
            Nonlinear_Solver* mMyNonLinSolver = nullptr;

            bool mIsMasterSystem = false;

            bool mBuildLinerSystemFlag = true;

            //! Map type. for special Petsc functionalities
            enum sol::MapType mMapType = sol::MapType::Epetra;

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
            Nonlinear_Problem( Solver_Interface* aSolverInterface,
                    const moris::sint            aNonlinearSolverManagerIndex = 0,
                    const bool                   aBuildLinerSystemFlag        = true,
                    const enum sol::MapType      aMapType                     = sol::MapType::Epetra );
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
            Nonlinear_Problem( sol::SOL_Warehouse* aNonlinDatabase,
                    Solver_Interface*              aSolverInterface,
                    sol::Dist_Vector*              aFullVector,
                    const moris::sint              aNonlinearSolverManagerIndex = 0,
                    const bool                     aBuildLinerSystemFlag        = true,
                    const enum sol::MapType        aMapType                     = sol::MapType::Epetra );

            //--------------------------------------------------------------------------------------------------
            ~Nonlinear_Problem();

            //--------------------------------------------------------------------------------------------------
            void set_interface( Solver_Interface* aSolverInterface );

            //--------------------------------------------------------------------------------------------------
            void
            set_nonlinear_solver( Nonlinear_Solver* aNonlinearSolver )
            {
                mMyNonLinSolver = aNonlinearSolver;
            };

            //--------------------------------------------------------------------------------------------------
            void build_linearized_problem( const bool& aRebuildJacobian,
                    const bool&                        aCombinedResJacAssebly,
                    sint                               aNonLinearIt );

            //--------------------------------------------------------------------------------------------------
            void build_linearized_problem( const bool& aRebuildJacobian,
                    const sint                         aNonLinearIt,
                    const sint                         aRestart );

            //--------------------------------------------------------------------------------------------------
            void print_sol_vec( const sint aNonLinearIt );

            //--------------------------------------------------------------------------------------------------
            void restart_from_sol_vec( const sint aNonLinearIt );

            //--------------------------------------------------------------------------------------------------
            dla::Linear_Problem*
            get_linearized_problem()
            {
                return mLinearProblem;
            };


            //--------------------------------------------------------------------------------------------------
            sol::Dist_Vector* get_full_vector();

            //--------------------------------------------------------------------------------------------------
            void extract_my_values( const moris::uint&      aNumIndices,
                    const moris::Matrix< DDSMat >&          aGlobalBlockRows,
                    const moris::uint&                      aBlockRowOffsets,
                    moris::Cell< moris::Matrix< DDRMat > >& LHSValues );

            //--------------------------------------------------------------------------------------------------
            void set_time_value( const moris::real& aLambda,
                    moris::uint                     aPos = 1 );
        };
    }    // namespace NLA
}    // namespace moris
#endif /* MORIS_DISTLINALG_CL_NLA_NONLINEAR_PROBLEM_HPP_ */

/*
 * cl_NLA_Nonlinear_Solver.hpp
 *
 *  Created on: Okt 6, 2018
 *      Author: schmidt
 */
#ifndef MORIS_DISTLINALG_CL_NLA_NONLINEAR_SOLVER_HPP_
#define MORIS_DISTLINALG_CL_NLA_NONLINEAR_SOLVER_HPP_

// MORIS header files.
#include "typedefs.hpp" // CON/src
#include "cl_Cell.hpp"
#include <memory>
#include "cl_Param_List.hpp"
#include "cl_MSI_Dof_Type_Enums.hpp"

#include "cl_NLA_Nonlinear_Solver_Enums.hpp"

namespace moris
{
    class Solver_Interface;
    namespace sol
    {
        class SOL_Warehouse;
    }
    namespace tsa
    {
        class Time_Solver_Algorithm;
    }
    namespace NLA
    {
        class Nonlinear_Problem;
        class Nonlinear_Algorithm;
        class Nonlinear_Solver
        {
            private:
                //! List of list of dof types
                moris::Cell< moris::Cell< enum MSI::Dof_Type > >  mStaggeredDofTypeList;

                //! List of secondary dependencies
                moris::Cell< moris::Cell< enum MSI::Dof_Type > >  mSecondaryDofTypeList;

                //! List with nonlinear solvers
                moris::Cell< std::shared_ptr< Nonlinear_Algorithm > > mNonlinearSolverAlgorithmList;

                //! List with nonlinear solvers
                moris::Cell< Nonlinear_Solver * > mNonLinearSubSolverList;

                //! Pointer to solver database
                sol::SOL_Warehouse * mSolverWarehouse = nullptr;

                //! Pointer to nonlinear problem
                Nonlinear_Problem * mNonlinearProblem = nullptr;

                //! Pointer to solver interface
                Solver_Interface * mSolverInput = nullptr;

                // Nonlinear solver manager index
                moris::sint mNonlinearSolverManagerIndex = -1;

                //! Reference norm
                moris::real mRefNorm = -1.0;

                //! First residual norm of this solve
                moris::real mFirstResidualNorm = -1.0;

                //! Actual residual norm
                moris::real mResidualNorm = -1.0;

                moris::ParameterList mParameterListNonLinearSolver;

                enum NonlinearSolverType mNonLinSolverType = NonlinearSolverType::END_ENUM;

                moris::uint mCallCounter = 0;
                moris::uint mCallCounterNonlinearSolver = 0;

                moris::sint mLevel = 0;

                moris::sint mTimeIter = 0;

                tsa::Time_Solver_Algorithm* mMyTimeSolverAlgorithm = nullptr;
                //        enum NonlinearSolverType mMyNonLinSolverType = NonlinearSolverType::END_ENUM;

            protected:

            public:
                /**
                 * @brief Constructor. Creates a default nonlinear solver
                 *
                 * @param[in] aNonLinSolverType Nonlinear solver type. Default is Newton
                 */
                Nonlinear_Solver( const enum NonlinearSolverType aNonLinSolverType = NonlinearSolverType::NEWTON_SOLVER );

                Nonlinear_Solver( const enum NonlinearSolverType aNonLinSolverType, const ParameterList aParameterlist );

                //--------------------------------------------------------------------------------------------------

                /**
                 * @brief Constructor. Sets given List of nonlinear solvers to this nonlinear solver manager
                 *
                 * @param[in] aNonlinerSolverList List of nonlinear solvers.
                 * @param[in] aNonLinSolverType Nonlinear solver type. Default is Newton
                 */
                Nonlinear_Solver(
                        moris::Cell< std::shared_ptr<Nonlinear_Algorithm > > & aNonlinerSolverList,
                        const enum NonlinearSolverType                         aNonLinSolverType = NonlinearSolverType::NEWTON_SOLVER );

                //--------------------------------------------------------------------------------------------------

                ~Nonlinear_Solver();

                //--------------------------------------------------------------------------------------------------

                void free_memory();

                //--------------------------------------------------------------------------------------------------

                /**
                 * @brief Sets one of the lists this nonlinear solver manager is operating on. Should be called multiple times for black solvers
                 *
                 * @param[in] aStaggeredDofTypeList List of dof types.
                 * @param[in] aLevel                Solver level in the block structure. Default is 0
                 */
                void set_dof_type_list(
                        const moris::Cell< enum MSI::Dof_Type > aStaggeredDofTypeList,
                        const moris::sint                       aLevel =  0);

                //--------------------------------------------------------------------------------------------------

                void set_secondary_dof_type_list( const moris::Cell< enum MSI::Dof_Type > aStaggeredDofTypeList);

                //--------------------------------------------------------------------------------------------------

                /**
                 * @brief Sets sub-non-linear solver this nonlinear solver is operating on
                 *
                 * @param[in] aNonlinearSolver Pointer to nonlinear solver
                 */
                void set_sub_nonlinear_solver( NLA::Nonlinear_Solver * aNonlinearSolver );

                //--------------------------------------------------------------------------------------------------

                /**
                 * @brief Sets sub-non-linear solver this nonlinear solver is operating on
                 *
                 * @param[in] aNonlinearSolver Pointer to nonlinear solver
                 * @param[in] aListEntry       List entry
                 */
                void set_sub_nonlinear_solver(
                        NLA::Nonlinear_Solver * aNonlinearSolver,
                        const moris::uint       aListEntry);

                //--------------------------------------------------------------------------------------------------

                /**
                 * @brief Sets solver interface
                 *
                 * @param[in] aSolverInterface Pointer to solver interface
                 */
                void set_solver_interface( Solver_Interface * aSolverInterface )
                {
                    mSolverInput = aSolverInterface;
                };

                //--------------------------------------------------------------------------------------------------

                /**
                 * @brief Gets solver interface
                 *
                 * @param[out] Pointer to solver interface
                 */
                Solver_Interface * get_solver_interface()
                {
                    return mSolverInput;
                };

                //--------------------------------------------------------------------------------------------------

                /**
                 * @brief Sets sub-non-linear solver this nonlinear solver is operating on
                 *
                 * @param[in] aNonlinearSolver Pointer to nonlinear solver
                 * @param[in] aListEntry       List entry
                 */
                Nonlinear_Solver * get_sub_nonlinear_solver( const moris::uint aListEntry)
                {
                    return mNonLinearSubSolverList( aListEntry );
                };

                //--------------------------------------------------------------------------------------------------

                /**
                 * @brief Set nonlinear solver. Uses push back to add the given nonlinear solver to the list.
                 *
                 * @param[in] aNonLinSolver Pointer to nonlinear solver.
                 */
                void set_nonlinear_algorithm( std::shared_ptr< Nonlinear_Algorithm > aNonLinSolver );

                //--------------------------------------------------------------------------------------------------

                /**
                 * @brief Set nonlinear solver on position in list
                 *
                 * @param[in] aNonLinSolver Pointer to nonlinear solver.
                 * @param[in] aListEntry Pointer to nonlinear solver.
                 */
                void set_nonlinear_algorithm(
                        std::shared_ptr< Nonlinear_Algorithm > aLinSolver,
                        const moris::uint                      aListEntry );

                //--------------------------------------------------------------------------------------------------

                /**
                 * @brief Get function for list of list of dof types
                 *
                 * @param[out] rListOfListsOfDofTypes Returns the nonlinear solver managers list of list of dof types
                 */
                moris::Cell< moris::Cell< enum MSI::Dof_Type > >  get_dof_type_list()
                {
                    return mStaggeredDofTypeList;
                };

                //--------------------------------------------------------------------------------------------------

                /**
                 * @brief Get the union of this nonlinear solver managers dof types
                 *
                 * @param[out] rUnionListOfDofTypes Returns the union list of this nonliner solver managers dof types
                 */
                moris::Cell< enum MSI::Dof_Type > get_dof_type_union();

                //--------------------------------------------------------------------------------------------------

                /**
                 * @brief Get secundary dof type list
                 *
                 * @param[out] rUnionListOfDofTypes Returns the union list of this nonliner solver managers dof types
                 */
                moris::Cell< moris::Cell< enum MSI::Dof_Type > >  get_secundary_dof_type_list()
                {
                    return mSecondaryDofTypeList;
                };

                //--------------------------------------------------------------------------------------------------

                moris::Cell< enum MSI::Dof_Type > get_sec_dof_type_union();

                //--------------------------------------------------------------------------------------------------

                /**
                 * @brief Sets the index of this nonlinear solver manager
                 *
                 * @param[in] aNonlinearSolverManagerIndex Nonlinear solver manager index
                 */
                void set_sonlinear_solver_manager_index( const moris::sint aNonlinearSolverManagerIndex );

                //--------------------------------------------------------------------------------------------------

                /**
                 * @brief Get the nonlinear solver manager index
                 *
                 * @param[out] rNonlinerSolverManagerIndex Returns the nonlinear solver manager index
                 */
                moris::sint get_sonlinear_solver_manager_index();

                //--------------------------------------------------------------------------------------------------

                /**
                 * @brief Retruns pointer to the solver database
                 *
                 * @param[out] rSolverDatabase Returns the pointer to the solver database
                 */
                sol::SOL_Warehouse * get_solver_warehouse(  )
                {
                    return mSolverWarehouse;
                };

                //--------------------------------------------------------------------------------------------------

                /**
                 * @brief Set pointer to the solver database
                 *
                 * @param[in] rSolverDatabase Poiner to the solver database
                 */
                void set_solver_warehouse( sol::SOL_Warehouse * aSolverWarehouse );

                //--------------------------------------------------------------------------------------------------
                /**
                 * @brief sets the time iteration
                 *
                 * @param[in] aTimeIter time iteration number
                 */
                void set_time_step_iter( const sint aTimeIter );

                //--------------------------------------------------------------------------------------------------
                /**
                 * @brief returns the time iteration value
                 */
                moris::sint get_time_step_iter(  );

                //--------------------------------------------------------------------------------------------------

                //void solve();

                void solve( sol::Dist_Vector * aFullVector);

                void solve( Nonlinear_Problem * aNonlinearProblem );

                void get_full_solution( moris::Matrix< DDRMat > & LHSValues );

                //--------------------------------------------------------------------------------------------------

                /**
                 * @brief Returns a reference to the reference norm
                 *
                 * @param[out] rRefNorm Returns the nonlinear solver manager reference norm
                 */
                moris::real & get_ref_norm(){ return mRefNorm; }

                //--------------------------------------------------------------------------------------------------

                /**
                 * @brief Returns a reference to the actual residual norm
                 *
                 * @param[out] rResidualNorm Returns the nonlinear solver manager residual norm
                 */
                moris::real & get_residual_norm(){ return mResidualNorm; }

                //--------------------------------------------------------------------------------------------------

                /**
                 * @brief Returns a reference to the actual first residual norm of this solve
                 *
                 * @param[out] rNonlinerSolverManagerIndex rResidualNorm Returns the first nonlinear solver manager residual norm
                 */
                moris::real & get_first_residual_norm()
                {
                    return mFirstResidualNorm;
                }

                //--------------------------------------------------------------------------------------------------

                void set_nonlinear_solver_manager_parameters();

                //--------------------------------------------------------------------------------------------------
                Nonlinear_Problem* get_my_nonlin_problem();

                void set_nonlin_solver_type( enum NonlinearSolverType aNonLinSolverType );

                enum NonlinearSolverType get_nonlin_solver_type();

                void set_time_solver_type( tsa::Time_Solver_Algorithm* aTimeSolverAlgorithm );

                //--------------------------------------------------------------------------------------------------

                ParameterListTypes & set_param( char const* aKey )
                {

                    return mParameterListNonLinearSolver( aKey );
                }
        };
    }
}
#endif /* MORIS_DISTLINALG_CL_NLA_NONLINEAR_SOLVER_HPP_ */

